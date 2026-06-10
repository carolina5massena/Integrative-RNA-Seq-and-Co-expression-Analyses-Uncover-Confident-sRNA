"""
pipeline_logic.py
-----------------
Parameterized implementations of the analysis steps, refactored so the PySide6
GUI can call them directly with user-customized inputs.

  - run_position_classification(params, log) -> step 01
  - run_wgcna(params, log)                   -> step 02 (calls Rscript)
  - run_deseq2(params, log)                  -> Differential Expression (Rscript)
  - run_filters(params, log)                 -> step 03

Each function accepts a `params` dict and a `log` callable (log(str)) used to
stream progress messages back to the GUI. They return the path of the main
output file produced.
"""

import glob
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd


# =====================================================================
# Rscript auto-detection (portable across machines)
# =====================================================================

def find_rscript():
    """Locate the Rscript executable without hardcoding a path.

    Resolution order:
      1. RSCRIPT environment variable (if it points to a real file)
      2. Rscript on the system PATH
      3. Common per-OS install locations (newest version wins on Windows)
      4. Fall back to the bare name "Rscript" (resolved at run time via PATH)
    """
    # 1. Explicit override
    env = os.environ.get("RSCRIPT", "").strip().strip('"')
    if env and os.path.isfile(env):
        return env

    # 2. On PATH
    found = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if found:
        return found

    # 3. Common install locations
    patterns = [
        # Windows
        r"C:\Program Files\R\*\bin\Rscript.exe",
        r"C:\Program Files\R\*\bin\x64\Rscript.exe",
        r"C:\Program Files (x86)\R\*\bin\Rscript.exe",
        os.path.expanduser(r"~\AppData\Local\Programs\R\*\bin\Rscript.exe"),
        # macOS
        "/Library/Frameworks/R.framework/Resources/bin/Rscript",
        "/opt/homebrew/bin/Rscript",
        "/usr/local/bin/Rscript",
        # Linux
        "/usr/bin/Rscript",
        "/usr/local/lib/R/bin/Rscript",
    ]
    candidates = []
    for pat in patterns:
        candidates.extend(glob.glob(pat))
    if candidates:
        # Prefer the highest version number (reverse sort puts R-4.6.0 > R-4.1.0)
        candidates.sort(reverse=True)
        return candidates[0]

    # 4. Last resort
    return "Rscript"


# =====================================================================
# 01 - POSITION CLASSIFICATION
# =====================================================================

def _read_gtf(gtf_path):
    cols = ["seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"]
    return pd.read_csv(gtf_path, sep="\t", comment="#", header=None, names=cols)


def _extract_attribute(attribute, key):
    if pd.isna(attribute):
        return None
    m = re.search(rf'{re.escape(key)} "([^"]+)"', attribute)
    return m.group(1) if m else None


def run_position_classification(p, log):
    """Classify sRNA coordinates relative to genes from a GTF file."""
    gtf_file = p["gtf_file"]
    log(f"Reading GTF: {gtf_file}")
    gtf = _read_gtf(gtf_file)
    gtf = gtf[gtf["feature"] == p["gtf_feature"]].copy()
    gtf["raw_id"] = gtf["attribute"].apply(lambda x: _extract_attribute(x, p["gtf_id_field"]))
    gtf = gtf[["raw_id", "start", "end", "strand"]].copy()
    gtf["length"] = gtf["end"] - gtf["start"]
    gtf = gtf[gtf["raw_id"].notna()].copy()

    srna_prefix = p["srna_prefix"]
    gene_prefix = p["gene_prefix"]
    filtered = gtf[gtf["raw_id"].astype(str).str.startswith((srna_prefix, gene_prefix))].copy()

    srn = filtered[filtered["raw_id"].str.startswith(srna_prefix)].copy().reset_index(drop=True)
    gene = filtered[filtered["raw_id"].str.startswith(gene_prefix)].copy().reset_index(drop=True)
    log(f"Found {len(srn)} sRNAs and {len(gene)} genes.")

    srn["exon_id"] = [f"{p['srna_rename_prefix']}{i+1}" for i in range(len(srn))]
    gene["gene_id"] = [f"{p['gene_rename_prefix']}{i+1}" for i in range(len(gene))]
    srn["original_id"] = srn["raw_id"]
    gene["original_id"] = gene["raw_id"]

    srn = srn[["exon_id", "original_id", "start", "end", "strand", "length"]].copy()
    gene = gene[["gene_id", "original_id", "start", "end", "strand", "length"]].copy()

    srn["location"] = None
    srn["gene_id"] = None
    srn[["start", "end"]] = srn[["start", "end"]].astype(int)
    gene[["start", "end"]] = gene[["start", "end"]].astype(int)

    overlap_thr = float(p["overlap_threshold"])
    utr = int(p["utr_distance"])

    log("Classifying Intragenic / Antisense...")
    for _, g in gene.iterrows():
        srn_len = srn["end"] - srn["start"]
        ov_start = srn["start"].clip(lower=g["start"])
        ov_end = srn["end"].clip(upper=g["end"])
        ov_len = (ov_end - ov_start).clip(lower=0)
        ratio = ov_len / srn_len
        mask = ratio >= overlap_thr
        same = srn["strand"] == g["strand"]
        srn.loc[mask & same, "location"] = "Intragenic"
        srn.loc[mask & ~same, "location"] = "Antisense"
        srn.loc[mask, "gene_id"] = g["gene_id"]

    log("Classifying UTRs...")
    for _, g in gene.iterrows():
        unclassified = srn["location"].isna()
        subset = srn[unclassified].copy()
        if subset.empty:
            break
        subset = subset[subset["strand"] == g["strand"]]
        if subset.empty:
            continue
        srn_len = subset["end"] - subset["start"]
        up_start = subset["start"].clip(lower=g["start"] - utr)
        up_end = subset["end"].clip(upper=g["start"])
        up_ratio = (up_end - up_start).clip(lower=0) / srn_len
        dn_start = subset["start"].clip(lower=g["end"])
        dn_end = subset["end"].clip(upper=g["end"] + utr)
        dn_ratio = (dn_end - dn_start).clip(lower=0) / srn_len
        up_mask = up_ratio >= overlap_thr
        dn_mask = dn_ratio >= overlap_thr
        if g["strand"] == "+":
            srn.loc[subset.index[up_mask], "location"] = "5' UTR"
            srn.loc[subset.index[dn_mask], "location"] = "3' UTR"
        else:
            srn.loc[subset.index[up_mask], "location"] = "3' UTR"
            srn.loc[subset.index[dn_mask], "location"] = "5' UTR"
        srn.loc[subset.index[up_mask | dn_mask], "gene_id"] = g["gene_id"]

    srn.loc[srn["location"].isna(), "location"] = "Intergenic"
    srn.loc[srn["location"] == "Intergenic", "gene_id"] = None

    out = p["output_file"]
    if out.lower().endswith(".csv"):
        srn.to_csv(out, index=False)
    else:
        srn.to_excel(out, index=False)
    log("Category counts:\n" + srn["location"].value_counts().to_string())
    log(f"Saved: {out}")
    return out


# =====================================================================
# Shared helpers
# =====================================================================

def _r_str(s):
    return '"' + str(s).replace("\\", "/").replace('"', '\\"') + '"'


def _safe_name(s):
    """Make a string safe to use inside a file name."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")
    return cleaned or "group"


def _run_rscript(script_path, rscript, cwd, log):
    log(f"Running: {rscript} {script_path}")
    proc = subprocess.Popen([rscript, script_path], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, cwd=cwd)
    for line in proc.stdout:
        log(line.rstrip())
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"Rscript exited with code {proc.returncode}")


# =====================================================================
# 02 - WGCNA CO-EXPRESSION
# =====================================================================

_WGCNA_TEMPLATE = r'''
library(WGCNA)
library(tibble)
library(fastDummies)
library(flashClust)

expr_counts_tsv   <- {expr_counts_tsv}
traits_csv        <- {traits_csv}
gene_id_col       <- {gene_id_col}
trait_sample_col  <- {trait_sample_col}
treatment_col     <- {treatment_col}
out_dir <- {out_dir}

out_sample_dendro_png        <- "sample_dendrogram_outlier_check.png"
out_sample_dendro_traits_png <- "sample_dendrogram_with_traits.png"
out_softpower_png            <- "soft_thresholding_power_selection.png"
out_gene_dendro_png          <- "gene_dendrogram_with_modules.png"
out_module_trait_pdf         <- "Module_Trait_Heatmap.pdf"
out_cyt_edges                <- "CytoscapeInput-edges.txt"
out_cyt_nodes                <- "CytoscapeInput-nodes.txt"
out_gene_module_csv          <- "gene_module_membership.csv"
out_top_hubs_csv             <- "top_hub_genes_per_module.csv"

n_threads <- {n_threads}
sample_tree_cut_height <- {sample_tree_cut_height}
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
softPower <- {softPower}
networkType <- {networkType}
adjacency_type <- {networkType}
corFnc <- {corFnc}
tomType <- {networkType}
minModuleSize <- {minModuleSize}
deepSplit <- {deepSplit}
pamRespectsDendro <- FALSE
MEDissThres <- {MEDissThres}
remove_first_dummy <- FALSE
remove_selected_columns <- TRUE
cytoscape_edge_threshold <- {cytoscape_edge_threshold}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
op <- function(x) file.path(out_dir, x)
enableWGCNAThreads(nThreads = n_threads)

raw_counts <- read.table(expr_counts_tsv, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
num_cols <- sapply(raw_counts, is.numeric)
num_cols[match(gene_id_col, names(raw_counts))] <- FALSE
SAbatch_df <- aggregate(raw_counts[, num_cols, drop = FALSE],
                        by = list(raw_counts[[gene_id_col]]), FUN = sum, na.rm = TRUE)
colnames(SAbatch_df)[1] <- gene_id_col
rownames(SAbatch_df) <- SAbatch_df[[gene_id_col]]
SAbatch <- SAbatch_df[, setdiff(names(SAbatch_df), gene_id_col), drop = FALSE]
SAdatExpr0 <- as.data.frame(t(SAbatch))

sampleTree <- hclust(dist(SAdatExpr0), method = "average")
png(op(out_sample_dendro_png), width = 12, height = 9, units = "in", res = 300)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = sample_tree_cut_height, col = "red")
dev.off()

datExpr <- SAdatExpr0
traitData <- read.csv(traits_csv, stringsAsFactors = FALSE)
Samples <- rownames(datExpr)
traitRows <- match(Samples, traitData[[trait_sample_col]])
datTraits <- traitData[traitRows, , drop = FALSE]
rownames(datTraits) <- datTraits[[trait_sample_col]]
datTraits_dummy <- fastDummies::dummy_cols(datTraits, select_columns = treatment_col,
                        remove_first_dummy = remove_first_dummy,
                        remove_selected_columns = remove_selected_columns)
datTraits_final <- datTraits_dummy[, !(names(datTraits_dummy) %in% trait_sample_col), drop = FALSE]
rownames(datTraits_final) <- rownames(datTraits)

sampleTree2 <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits_final, signed = FALSE)
png(op(out_sample_dendro_traits_png), width = 12, height = 9, units = "in", res = 300)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits_final),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

sft <- pickSoftThreshold(datExpr, networkType = networkType, powerVector = powers, verbose = 5)
png(op(out_softpower_png), width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2)); cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

adjacency <- adjacency(datExpr, power = softPower, type = adjacency_type, corFnc = corFnc)
TOM <- TOMsimilarity(adjacency, TOMType = tomType)
dissTOM <- 1 - TOM
geneTree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = deepSplit,
                             pamRespectsDendro = pamRespectsDendro, minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors <- merge$colors
MEs <- merge$newMEs

png(op(out_gene_dendro_png), width = 12, height = 9, units = "in", res = 300)
plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors),
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

nSamples <- nrow(datExpr)
moduleTraitCor <- cor(MEs, datTraits_final, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
pdf(file = op(out_module_trait_pdf), width = 10, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits_final), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5,
               zlim = c(-1, 1), main = "Module-trait relationships")
dev.off()

probes <- names(datExpr)
inModule <- is.finite(match(moduleColors, unique(moduleColors[moduleColors != "grey"])))
modProbes <- probes[inModule]
modTOM <- dissTOM[inModule, inModule]
diag(modTOM) <- 0
exportNetworkToCytoscape(1 - modTOM, edgeFile = op(out_cyt_edges), nodeFile = op(out_cyt_nodes),
                         weighted = TRUE, threshold = cytoscape_edge_threshold,
                         nodeNames = modProbes, nodeAttr = moduleColors[inModule])
geneModuleMembership <- data.frame(Gene = probes, Module = moduleColors)
write.csv(geneModuleMembership, op(out_gene_module_csv), row.names = FALSE)
topHubs <- chooseTopHubInEachModule(datExpr, colorh = moduleColors, omitColors = "grey", power = softPower)
write.csv(topHubs, op(out_top_hubs_csv), row.names = FALSE)
message("WGCNA analysis complete. Check outputs in: ", normalizePath(out_dir))
'''


def run_wgcna(p, log):
    """Render the WGCNA R script with user params and run it via Rscript."""
    nthreads = str(p.get("n_threads", "")).strip()
    n_threads_r = nthreads if nthreads else "NULL"
    script = _WGCNA_TEMPLATE.format(
        expr_counts_tsv=_r_str(p["expr_counts_tsv"]),
        traits_csv=_r_str(p["traits_csv"]),
        gene_id_col=_r_str(p["gene_id_col"]),
        trait_sample_col=_r_str(p["trait_sample_col"]),
        treatment_col=_r_str(p["treatment_col"]),
        out_dir=_r_str(p["out_dir"]),
        n_threads=n_threads_r,
        sample_tree_cut_height=p["sample_tree_cut_height"],
        softPower=p["softPower"],
        networkType=_r_str(p["networkType"]),
        corFnc=_r_str(p["corFnc"]),
        minModuleSize=p["minModuleSize"],
        deepSplit=p["deepSplit"],
        MEDissThres=p["MEDissThres"],
        cytoscape_edge_threshold=p["cytoscape_edge_threshold"],
    )
    script_path = os.path.join(p["out_dir"], "WGCNA_configured.R")
    with open(script_path, "w") as f:
        f.write(script)
    log(f"Wrote configured R script: {script_path}")
    rscript = str(p.get("rscript_path", "")).strip() or "Rscript"
    _run_rscript(script_path, rscript, p["out_dir"], log)
    return os.path.join(p["out_dir"], "CytoscapeInput-nodes.txt")


# =====================================================================
# DIFFERENTIAL EXPRESSION (DESeq2 in R)
# =====================================================================

_DESEQ2_TEMPLATE = r'''
suppressMessages(library(DESeq2))

counts_tsv   <- {counts_tsv}
meta_csv     <- {meta_csv}
gene_id_col  <- {gene_id_col}
sample_col   <- {sample_col}
condition_col<- {condition_col}
ref_level    <- {ref_level}
test_level   <- {test_level}
min_count    <- {min_count}
out_tsv      <- {out_tsv}

counts <- read.table(counts_tsv, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- counts[[gene_id_col]]
counts[[gene_id_col]] <- NULL
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

coldata <- read.csv(meta_csv, stringsAsFactors = FALSE, check.names = FALSE)
rownames(coldata) <- coldata[[sample_col]]

common <- intersect(colnames(counts), rownames(coldata))
if (length(common) < 2) stop("Fewer than 2 matching samples between counts and metadata.")
counts  <- counts[, common, drop = FALSE]
coldata <- coldata[common, , drop = FALSE]

coldata[[condition_col]] <- factor(coldata[[condition_col]])
coldata[[condition_col]] <- relevel(coldata[[condition_col]], ref = ref_level)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = as.formula(paste0("~", condition_col)))
dds <- dds[rowSums(counts(dds)) >= min_count, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c(condition_col, test_level, ref_level))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", "log2FoldChange", "padj",
                     "baseMean", "pvalue", "lfcSE", "stat")]

write.table(res_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("DESeq2 complete. ", nrow(res_df), " genes written to ", out_tsv)
message("Significant (padj < 0.05): ", sum(res_df$padj < 0.05, na.rm = TRUE))
'''


# Per-comparison template. Estimates a shared dispersion from ALL samples
# (the strains act as replicates of each condition), then runs one Wald test
# per comparison group, writing <out_base>_<group><out_ext> for each. Groups
# that happen to have their own replicates use their own dispersion instead.
_DESEQ2_COMPARISON_TEMPLATE = r'''
suppressMessages(library(DESeq2))

counts_tsv     <- {counts_tsv}
meta_csv       <- {meta_csv}
gene_id_col    <- {gene_id_col}
sample_col     <- {sample_col}
condition_col  <- {condition_col}
comparison_col <- {comparison_col}
ref_level      <- {ref_level}
test_level     <- {test_level}
min_count      <- {min_count}
out_dir        <- {out_dir}
out_base       <- {out_base}
out_ext        <- {out_ext}

counts <- read.table(counts_tsv, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- counts[[gene_id_col]]
counts[[gene_id_col]] <- NULL
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

coldata <- read.csv(meta_csv, stringsAsFactors = FALSE, check.names = FALSE)
rownames(coldata) <- coldata[[sample_col]]

common <- intersect(colnames(counts), rownames(coldata))
if (length(common) < 2) stop("Fewer than 2 matching samples between counts and metadata.")
counts  <- counts[, common, drop = FALSE]
coldata <- coldata[common, , drop = FALSE]

coldata[[condition_col]] <- relevel(factor(coldata[[condition_col]]), ref = ref_level)

## ---- shared dispersion estimated from ALL samples (strains as replicates) ----
dds_full <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                   design = as.formula(paste0("~", condition_col)))
dds_full <- dds_full[rowSums(counts(dds_full)) >= min_count, ]
dds_full <- estimateSizeFactors(dds_full)
dds_full <- estimateDispersions(dds_full)
gene_universe <- rownames(dds_full)
shared_disp <- dispersions(dds_full)
names(shared_disp) <- gene_universe
message("Shared dispersion estimated from ", ncol(dds_full), " samples over ",
        length(gene_universe), " genes.")

cols <- c("gene_id", "log2FoldChange", "padj", "baseMean", "pvalue", "lfcSE", "stat")
groups <- unique(as.character(coldata[[comparison_col]]))
groups <- groups[!is.na(groups)]

for (g in groups) {{
  cd <- coldata[as.character(coldata[[comparison_col]]) == g, , drop = FALSE]
  conds <- as.character(cd[[condition_col]])
  if (!(ref_level %in% conds) || !(test_level %in% conds)) {{
    message("[skip] comparison ", g, ": needs both ", ref_level, " and ", test_level)
    next
  }}
  cnt <- counts[gene_universe, rownames(cd), drop = FALSE]
  cd[[condition_col]] <- relevel(factor(cd[[condition_col]]), ref = ref_level)
  dds <- DESeqDataSetFromMatrix(countData = cnt, colData = cd,
                                design = as.formula(paste0("~", condition_col)))
  dds <- estimateSizeFactors(dds)
  if (min(table(cd[[condition_col]])) >= 2) {{
    dds <- estimateDispersions(dds)
    note <- "own dispersion"
  }} else {{
    dispersionFunction(dds) <- dispersionFunction(dds_full)
    dispersions(dds) <- shared_disp[gene_universe]
    note <- "shared dispersion (no replicates)"
  }}
  dds <- nbinomWaldTest(dds)
  res <- results(dds, contrast = c(condition_col, test_level, ref_level))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, cols]
  safe_g <- gsub("[^A-Za-z0-9._-]+", "_", g)
  out_tsv <- file.path(out_dir, paste0(out_base, "_", safe_g, out_ext))
  write.table(res_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  message("comparison ", g, " -> ", basename(out_tsv), " [", note, "]  sig padj<0.05: ",
          sum(res_df$padj < 0.05, na.rm = TRUE))
}}
message("All comparisons complete.")
'''


def _run_one_deseq2(p, meta_csv, out_tsv, out_dir, rscript, log):
    """Render the DESeq2 R script for one metadata table and run it."""
    script = _DESEQ2_TEMPLATE.format(
        counts_tsv=_r_str(p["counts_tsv"]),
        meta_csv=_r_str(meta_csv),
        gene_id_col=_r_str(p["gene_id_col"]),
        sample_col=_r_str(p["sample_col"]),
        condition_col=_r_str(p["condition_col"]),
        ref_level=_r_str(p["ref_level"]),
        test_level=_r_str(p["test_level"]),
        min_count=int(p["min_count"]),
        out_tsv=_r_str(out_tsv),
    )
    tag = _safe_name(os.path.splitext(os.path.basename(out_tsv))[0])
    script_path = os.path.join(out_dir, f"DESeq2_configured_{tag}.R")
    with open(script_path, "w") as f:
        f.write(script)
    log(f"Wrote configured R script: {script_path}")
    _run_rscript(script_path, rscript, out_dir, log)


def run_deseq2(p, log):
    """Render a DESeq2 R script with user params and run it via Rscript.

    If a `comparison_col` is given, the metadata is split by that column and a
    separate DESeq2 run (reference vs test) is performed for each group,
    writing one DEG TSV per group (e.g. DEG_result_A.tsv, DEG_result_B.tsv).
    The returned value is a comma-separated list of every output file, ready to
    paste straight into the 03 Filters tab for a consensus.
    """
    out_tsv = p["output_tsv"]
    out_dir = os.path.dirname(out_tsv) or "."
    rscript = str(p.get("rscript_path", "")).strip() or "Rscript"
    comparison_col = str(p.get("comparison_col", "")).strip()
    sample_col = p["sample_col"]
    condition_col = p["condition_col"]
    ref_level = str(p["ref_level"])
    test_level = str(p["test_level"])

    # --- Single combined run (original behaviour) ---
    if not comparison_col:
        _run_one_deseq2(p, p["meta_csv"], out_tsv, out_dir, rscript, log)
        return out_tsv

    # --- One run per comparison group ---
    meta = pd.read_csv(p["meta_csv"])
    if comparison_col not in meta.columns:
        raise ValueError(
            f"Comparison column '{comparison_col}' not found in metadata. "
            f"Available columns: {list(meta.columns)}"
        )
    groups = list(pd.unique(meta[comparison_col].dropna()))
    if not groups:
        raise ValueError(f"No values found in comparison column '{comparison_col}'.")
    log(f"Comparison column '{comparison_col}' -> {len(groups)} group(s): "
        f"{', '.join(map(str, groups))}")

    base, ext = os.path.splitext(out_tsv)
    ext = ext or ".tsv"

    # Validate groups up front (the R script also skips invalid ones).
    valid = []
    for g in groups:
        sub = meta[meta[comparison_col] == g]
        conds = set(sub[condition_col].astype(str)) if condition_col in sub.columns else set()
        if ref_level not in conds or test_level not in conds:
            log(f"  [skip] comparison '{g}': needs both '{ref_level}' and "
                f"'{test_level}' in '{condition_col}', found {sorted(conds)}")
            continue
        valid.append(g)

    if not valid:
        raise RuntimeError(
            "No comparison groups could be run. Check that each group in "
            f"'{comparison_col}' contains both the reference ('{ref_level}') and "
            f"test ('{test_level}') levels of '{condition_col}'."
        )

    log("Note: comparisons without within-group replicates use a shared "
        "dispersion estimated from all samples (DESeq2 cannot estimate "
        "dispersion from a single sample per condition). Groups that do have "
        "replicates use their own dispersion.")

    script = _DESEQ2_COMPARISON_TEMPLATE.format(
        counts_tsv=_r_str(p["counts_tsv"]),
        meta_csv=_r_str(p["meta_csv"]),
        gene_id_col=_r_str(p["gene_id_col"]),
        sample_col=_r_str(sample_col),
        condition_col=_r_str(condition_col),
        comparison_col=_r_str(comparison_col),
        ref_level=_r_str(ref_level),
        test_level=_r_str(test_level),
        min_count=int(p["min_count"]),
        out_dir=_r_str(out_dir),
        out_base=_r_str(os.path.basename(base)),
        out_ext=_r_str(ext),
    )
    script_path = os.path.join(out_dir, "DESeq2_configured_comparisons.R")
    with open(script_path, "w") as f:
        f.write(script)
    log(f"Wrote configured R script: {script_path}")
    _run_rscript(script_path, rscript, out_dir, log)

    outputs = [f"{base}_{_safe_name(g)}{ext}" for g in valid]
    log(f"Done. {len(outputs)} DEG file(s) written:")
    for o in outputs:
        log(f"  {o}")
    log("Tip: paste all of these into 'DEG TSV file(s)' on the 03 Filters tab "
        "to run a consensus across comparisons.")
    return ", ".join(outputs)


# =====================================================================
# 03 - FILTERS  (fully customizable)
# =====================================================================

def _check_exists(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")


# Definition of every energy/probability metric:
#   key -> (direction, default_threshold)
# direction "max": keep rows where value <= threshold
# direction "min": keep rows where value >= threshold
_ENERGY_METRICS = {
    "E_intaRNA": ("max", -2.44),
    "E_Rnaplex": ("max", -32.6),
    "E_TargetRNA3": ("max", -5.13),
    "Probability_TargetRNA3": ("min", 0.06),
    "Probability_sRNARFTarget": ("min", 0.40),
}


def _build_deg_status(deg_files, padj_cutoff, min_strains, log, min_basemean=0.0):
    """Return a dict gene_id -> 'upregulated'/'downregulated' from DEG TSVs.

    A gene passes within a file when padj <= padj_cutoff (and, if min_basemean > 0
    and the file has a baseMean column, baseMean >= min_basemean).
    """
    deg_tables = {}
    for fp in deg_files:
        _check_exists(fp)
        df = pd.read_csv(fp, sep="\t")
        if not {"gene_id", "log2FoldChange", "padj"}.issubset(df.columns):
            raise ValueError(f"{fp} must contain gene_id, log2FoldChange, padj")
        sel = df["padj"] <= padj_cutoff
        if min_basemean and min_basemean > 0:
            if "baseMean" in df.columns:
                sel &= df["baseMean"] >= min_basemean
            else:
                log(f"  [warn] {os.path.basename(fp)} has no baseMean column; "
                    f"baseMean filter skipped for it.")
        deg_tables[fp] = df.loc[sel, ["gene_id", "log2FoldChange"]].copy()
    merged = None
    for i, (fp, df) in enumerate(deg_tables.items()):
        tmp = df.rename(columns={"log2FoldChange": f"strain{i}"})
        merged = tmp if merged is None else pd.merge(merged, tmp, on="gene_id", how="outer")
    if merged is None:
        return {}
    num = merged.select_dtypes(include=[np.number])

    def row_status(r):
        if (r > 0).sum() >= min_strains:
            return "upregulated"
        if (r < 0).sum() >= min_strains:
            return "downregulated"
        return np.nan

    merged["regulation_status"] = num.apply(row_status, axis=1)
    deg_cons = merged.dropna(subset=["regulation_status"])
    log(f"Consistent DEGs: {len(deg_cons)}")
    return deg_cons.set_index("gene_id")["regulation_status"].to_dict()


def run_filters(p, log):
    """Energy/probability + DEG + WGCNA-module + network-weight filtering.

    All filters are optional/configurable via the params dict (sensible defaults
    keep the original behaviour). Recognized keys:

      predictions_csv, deg_files, module_nodes_file, module_edges_file, output_csv
                             (any empty path skips that step; if predictions_csv is
                              empty the pairs are built from module_edges_file)
      srna_prefix            : str   (node-name prefix that marks an sRNA when
                                      building pairs from edges; default 'sRNA')
      filters_outlier        : {metric: threshold}
      energy_enabled         : {metric: bool}            (default all True)
                             (a metric is also skipped automatically if its column
                              is absent from the predictions file)
      energy_min_pass        : int  (how many enabled energy filters must pass;
                                     default = number of enabled = logical AND)
      padj_cutoff            : float
      min_strains_consistent : int
      deg_min_basemean       : float (min DESeq2 baseMean per DEG; 0 = off)
      require_srna_deg       : bool   (sRNA must be a consistent DEG)
      require_target_deg     : bool   (target must be a consistent DEG)
      srna_direction         : 'any'|'upregulated'|'downregulated'
      target_direction       : 'any'|'upregulated'|'downregulated'
      pair_relationship      : 'any'|'same'|'opposite'   (sRNA vs target direction)
      module_mode            : 'same'|'different'|'any'
      require_in_network      : bool   (interaction must have a WGCNA edge weight)
      min_edge_weight        : float
      selection_mode         : 'best_per_target'|'best_per_srna'|'all'
    """
    # ---- 0. Base table of sRNA-Target pairs --------------------------------
    # Prediction CSV is preferred (it carries the program scores). If it is left
    # empty, the pairs are built from the WGCNA edges file instead, and the
    # energy/probability step is skipped (no scores available).
    pred_path = str(p.get("predictions_csv", "")).strip()
    edges_path = str(p.get("module_edges_file", "")).strip()

    if pred_path:
        _check_exists(pred_path)
        pred = pd.read_csv(pred_path)
        if not {"sRNA", "Target"}.issubset(pred.columns):
            raise ValueError("Prediction file must contain at least 'sRNA' and 'Target' columns.")
        log(f"Total predictions: {len(pred)}")
        have_scores = True
    else:
        if not edges_path:
            raise ValueError(
                "No Prediction CSV and no WGCNA edges file were given. Provide a "
                "Prediction CSV (scored sRNA-target pairs) or an edges file to "
                "build the pairs from the co-expression network."
            )
        _check_exists(edges_path)
        srna_prefix = str(p.get("srna_prefix", "sRNA")).strip()
        edf = pd.read_csv(edges_path, sep="\t")
        if not {"fromNode", "toNode"}.issubset(edf.columns):
            raise ValueError("edges file must contain fromNode and toNode columns.")
        from_s = edf["fromNode"].astype(str).str.startswith(srna_prefix)
        to_s = edf["toNode"].astype(str).str.startswith(srna_prefix)
        keep = from_s ^ to_s  # exactly one endpoint is an sRNA
        edf = edf.loc[keep].copy()
        if edf.empty:
            raise ValueError(
                f"No sRNA-target edges found: no edges have exactly one node whose "
                f"name starts with the sRNA prefix '{srna_prefix}'. Adjust the prefix."
            )
        edf["sRNA"] = np.where(from_s.loc[keep], edf["fromNode"], edf["toNode"])
        edf["Target"] = np.where(from_s.loc[keep], edf["toNode"], edf["fromNode"])
        pred = edf[["sRNA", "Target"]].copy()
        log(f"No Prediction CSV given - built {len(pred)} sRNA-target pair(s) from the "
            f"edges file (sRNA prefix '{srna_prefix}').")
        have_scores = False

    # ---- 1. Energy / probability filters (optional; skips missing columns) --
    if have_scores:
        thr = p.get("filters_outlier", {})
        enabled = p.get("energy_enabled", {m: True for m in _ENERGY_METRICS})
        pass_cols = []
        used = []
        for metric, (direction, default_thr) in _ENERGY_METRICS.items():
            if not enabled.get(metric, True):
                continue
            if metric not in pred.columns:
                log(f"  [skip] {metric}: not present in the Prediction file.")
                continue
            t = thr.get(metric + ("_max" if direction == "max" else "_min"),
                        thr.get(metric, default_thr))
            if direction == "max":
                pass_cols.append(pred[metric] <= t)
            else:
                pass_cols.append(pred[metric] >= t)
            used.append(metric)
        if pass_cols:
            passes = pd.concat(pass_cols, axis=1).sum(axis=1)
            n_enabled = len(pass_cols)
            min_pass = int(p.get("energy_min_pass", n_enabled) or n_enabled)
            min_pass = max(1, min(min_pass, n_enabled))
            pred = pred.loc[passes >= min_pass].copy()
            log(f"After energy/probability filter (>= {min_pass}/{n_enabled} of {used}): {len(pred)}")
        else:
            log("Energy/probability filter skipped (no enabled metric columns present).")
    else:
        log("Energy/probability filter skipped (pairs built from edges; no scores).")

    # ---- 2. DEG status + direction -----------------------------------------
    deg_map = _build_deg_status(p["deg_files"], p.get("padj_cutoff", 0.05),
                                p.get("min_strains_consistent", 1), log,
                                p.get("deg_min_basemean", 0.0)) if p.get("deg_files") else {}
    pred["sRNA_DEG"] = pred["sRNA"].map(deg_map)
    pred["Target_DEG"] = pred["Target"].map(deg_map)

    if p.get("require_srna_deg", True):
        pred = pred.dropna(subset=["sRNA_DEG"])
    if p.get("require_target_deg", True):
        pred = pred.dropna(subset=["Target_DEG"])

    sdir = p.get("srna_direction", "any")
    if sdir != "any":
        pred = pred[pred["sRNA_DEG"] == sdir]
    tdir = p.get("target_direction", "any")
    if tdir != "any":
        pred = pred[pred["Target_DEG"] == tdir]

    rel = p.get("pair_relationship", "any")
    if rel == "same":
        pred = pred[pred["sRNA_DEG"] == pred["Target_DEG"]]
    elif rel == "opposite":
        pred = pred[pred["sRNA_DEG"].notna() & pred["Target_DEG"].notna()
                    & (pred["sRNA_DEG"] != pred["Target_DEG"])]
    log(f"After DEG filter: {len(pred)}")

    # ---- 3. WGCNA module relationship --------------------------------------
    module_mode = p.get("module_mode", "same")
    if module_mode != "any" and p.get("module_nodes_file"):
        _check_exists(p["module_nodes_file"])
        nodes = pd.read_csv(p["module_nodes_file"], sep="\t")
        nodes = nodes.rename(columns={"nodeAttr[nodesPresent, ]": "nodeAttr"})
        if not {"nodeName", "nodeAttr"}.issubset(nodes.columns):
            raise ValueError("nodes file must contain nodeName, nodeAttr")
        nm = nodes.set_index("nodeName")["nodeAttr"].to_dict()
        pred["sRNA_Module"] = pred["sRNA"].map(nm)
        pred["Target_Module"] = pred["Target"].map(nm)
        both = pred["sRNA_Module"].notna() & pred["Target_Module"].notna()
        if module_mode == "same":
            pred = pred[both & (pred["sRNA_Module"] == pred["Target_Module"])]
        elif module_mode == "different":
            pred = pred[both & (pred["sRNA_Module"] != pred["Target_Module"])]
        log(f"After module filter ({module_mode}): {len(pred)}")
    else:
        log("Module filter skipped.")

    # ---- 4. Network edge weight + selection ---------------------------------
    if p.get("module_edges_file"):
        _check_exists(p["module_edges_file"])
        edges = pd.read_csv(p["module_edges_file"], sep="\t")
        if not {"fromNode", "toNode", "weight"}.issubset(edges.columns):
            raise ValueError("edges file must contain fromNode, toNode, weight")
        edges["pair"] = edges.apply(lambda r: tuple(sorted([r["fromNode"], r["toNode"]])), axis=1)
        edges = edges[["pair", "weight"]]
        pred = pred.copy()
        pred["pair"] = pred.apply(lambda r: tuple(sorted([r["sRNA"], r["Target"]])), axis=1)
        pred = pred.merge(edges, on="pair", how="left").drop(columns=["pair"])
    else:
        pred = pred.copy()
        pred["weight"] = np.nan

    if p.get("require_in_network", True):
        pred = pred[pred["weight"].notna()]
        log(f"After requiring network edge: {len(pred)}")
    min_w = float(p.get("min_edge_weight", 0.0) or 0.0)
    if min_w > 0:
        pred = pred[(pred["weight"].isna()) | (pred["weight"] >= min_w)]
        log(f"After min edge weight ({min_w}): {len(pred)}")

    selection = p.get("selection_mode", "best_per_target")
    if selection in ("best_per_target", "best_per_srna"):
        key = "Target" if selection == "best_per_target" else "sRNA"
        pred = pred.sort_values([key, "weight"], ascending=[True, False])
        pred = pred.drop_duplicates(key)
        log(f"After {selection} selection: {len(pred)}")

    out = p["output_csv"]
    pred.to_csv(out, index=False)
    log(f"Final interactions: {len(pred)}")
    log(f"Saved: {out}")
    return out
