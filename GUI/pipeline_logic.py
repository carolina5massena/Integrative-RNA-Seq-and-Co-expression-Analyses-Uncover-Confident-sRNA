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
    paste straight into the 05 Filters tab for a consensus.
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
    log("Tip: paste all of these into 'DEG TSV file(s)' on the 05 Filters tab "
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
    "p_intaRNA": ("max", 0.41),
    "p_Rnaplex": ("max", 0.76),
    "p_TargetRNA3": ("max", 0.89),
}


def _build_deg_status(deg_files, padj_cutoff, min_strains, log,
                      min_basemean=0.0, basemean_enabled=False,
                      lfc_up_enabled=False, lfc_up=1.0,
                      lfc_down_enabled=False, lfc_down=-1.0,
                      consensus_mode="sign"):
    """Build DEG gene sets and a consensus status map from the DEG TSVs.

    A gene is "called" within one file when it passes, in this cascade:
        1. padj     <= padj_cutoff
        2. baseMean >= min_basemean          (only if basemean_enabled and the file
                                              has a baseMean column)
        3. log2FC outside the inclusion band (only if a log2FC side is enabled):
               log2FC >= lfc_up  OR  log2FC <= lfc_down
    A gene absent from a file is simply not called there (as if padj were infinite,
    baseMean negative, log2FC 0 for that file).

    Consensus needs the call in at least `min_strains` (X) files:
        consensus_mode == "sign"    : up   if log2FC > 0 in >= X files,
                                      down if log2FC < 0 in >= X files
        consensus_mode == "presence": kept if the gene is called in >= X files at
                                      all (log2FC sign ignored); direction is then
                                      set by the majority sign, for the optional
                                      direction sub-filters.

    Returns a dict:
        genes_padj      : set  (padj passed in >= 1 file)
        genes_basemean  : set  (padj + baseMean passed in >= 1 file)
        genes_log2fc    : set  (padj + baseMean + log2FC band passed in >= 1 file)
        status_map      : {gene: 'upregulated'|'downregulated'}  (consensus genes)
        log2fc_map      : {gene: mean log2FoldChange across the files that called it}
    """
    from collections import defaultdict
    padj_count = defaultdict(int)   # files where padj passed
    bm_count = defaultdict(int)     # files where padj + baseMean passed
    call_count = defaultdict(int)   # files where the full per-file DEG call passed
    pos_count = defaultdict(int)    # called files with log2FC > 0
    neg_count = defaultdict(int)    # called files with log2FC < 0
    lfc_sum = defaultdict(float)    # sum of log2FC over called files

    band_on = bool(lfc_up_enabled) or bool(lfc_down_enabled)

    for fp in deg_files:
        _check_exists(fp)
        df = pd.read_csv(fp, sep="\t")
        if not {"gene_id", "log2FoldChange", "padj"}.issubset(df.columns):
            raise ValueError(f"{fp} must contain gene_id, log2FoldChange, padj")
        # Tolerate text / comma-decimal numbers in the numeric columns.
        for c in ("padj", "log2FoldChange", "baseMean"):
            if c in df.columns:
                df[c] = _numeric_col(df[c])

        padj_ok = df["padj"] <= padj_cutoff
        bm_ok = padj_ok.copy()
        if basemean_enabled and min_basemean and min_basemean > 0:
            if "baseMean" in df.columns:
                bm_ok = bm_ok & (df["baseMean"] >= min_basemean)
            else:
                log(f"  [warn] {os.path.basename(fp)} has no baseMean column; "
                    f"baseMean filter skipped for it.")
        if band_on:
            band = pd.Series(False, index=df.index)
            if lfc_up_enabled:
                band = band | (df["log2FoldChange"] >= lfc_up)
            if lfc_down_enabled:
                band = band | (df["log2FoldChange"] <= lfc_down)
            call_ok = bm_ok & band
        else:
            call_ok = bm_ok

        for g in df.loc[padj_ok, "gene_id"]:
            padj_count[g] += 1
        for g in df.loc[bm_ok, "gene_id"]:
            bm_count[g] += 1
        for g, l in zip(df.loc[call_ok, "gene_id"], df.loc[call_ok, "log2FoldChange"]):
            call_count[g] += 1
            lfc_sum[g] += float(l)
            if l > 0:
                pos_count[g] += 1
            elif l < 0:
                neg_count[g] += 1

    genes_padj = set(padj_count)
    genes_basemean = set(bm_count)
    genes_log2fc = set(call_count)

    X = int(min_strains)
    status_map = {}
    log2fc_map = {}
    for g, c in call_count.items():
        if consensus_mode == "presence":
            if c >= X:
                status_map[g] = ("upregulated" if pos_count[g] >= neg_count[g]
                                 else "downregulated")
        else:  # "sign"
            if pos_count[g] >= X:
                status_map[g] = "upregulated"
            elif neg_count[g] >= X:
                status_map[g] = "downregulated"
        if g in status_map:
            log2fc_map[g] = lfc_sum[g] / c if c else 0.0

    log(f"  DEG gene sets: padj={len(genes_padj)}, +baseMean={len(genes_basemean)}, "
        f"+log2FC band={len(genes_log2fc)}, "
        f"consensus[{consensus_mode}, >={X} file(s)]={len(status_map)}")
    return {
        "genes_padj": genes_padj,
        "genes_basemean": genes_basemean,
        "genes_log2fc": genes_log2fc,
        "status_map": status_map,
        "log2fc_map": log2fc_map,
    }


def _both_in(pred, gene_set):
    """Boolean mask: pairs whose sRNA AND target are both in gene_set."""
    return pred["sRNA"].isin(gene_set) & pred["Target"].isin(gene_set)


def _numeric_col(s):
    """Coerce a column to numeric, tolerating text and comma decimals.

    Values that are already numeric pass straight through. Text values (as
    produced by a pt-BR / Excel export, e.g. "-69,9" or "3,0e-05", or stray
    entries like "NA") are parsed; anything that still cannot be understood
    becomes NaN. This prevents a "'<=' not supported between str and float"
    crash when a score column is read as text.
    """
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")
    txt = s.astype(str).str.strip()
    out = pd.to_numeric(txt, errors="coerce")
    # Retry the ones that failed, reading comma as the decimal separator and
    # dropping any thousands dots ("1.234,5" -> "1234.5", "-69,9" -> "-69.9").
    mask = out.isna() & txt.ne("") & ~txt.str.lower().isin(["nan", "na", "none"])
    if mask.any():
        fixed = (txt[mask].str.replace(".", "", regex=False)
                          .str.replace(",", ".", regex=False))
        out.loc[mask] = pd.to_numeric(fixed, errors="coerce")
    return out


def _pair_key(a, b):
    """Order-independent string key for a node pair, vectorized over two Series.

    Returns "lo\\thi" so (X, Y) and (Y, X) map to the same key. Much faster and
    lighter than a row-wise apply(sorted(...)) when the edges file has millions of
    rows.
    """
    a = a.astype(str).to_numpy()
    b = b.astype(str).to_numpy()
    swap = a > b
    lo = np.where(swap, b, a)
    hi = np.where(swap, a, b)
    return pd.Series([f"{x}\t{y}" for x, y in zip(lo, hi)])


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
      lfc_up_enabled         : bool   (enable the 'log2FC >= lfc_up' side of the band)
      lfc_up                 : float  (keep genes with log2FC >= this)
      lfc_down_enabled       : bool   (enable the 'log2FC <= lfc_down' side)
      lfc_down               : float  (keep genes with log2FC <= this)
                             (band applies to BOTH sRNA and target; a gene passes
                              when log2FC >= lfc_up OR <= lfc_down. Both disabled =
                              no DEG requirement on the endpoints.)
      srna_direction         : 'any'|'upregulated'|'downregulated'
      target_direction       : 'any'|'upregulated'|'downregulated'
      pair_relationship      : 'any'|'same'|'opposite'   (sRNA vs target direction)
      location_filter_enabled: bool   (apply the sRNA genomic location filter)
      srna_annotation_file   : str    (annotation table from tab 01; .xlsx/.csv)
      location_id_col        : str    (sRNA ID column in the annotation; 'exon_id')
      locations_keep         : list   (categories to keep, e.g. ['Intragenic',
                                       'Antisense', "5' UTR", "3' UTR", 'Intergenic'])
      module_mode            : 'same'|'different'|'any'
      require_in_network      : bool   (interaction must have a WGCNA edge weight)
      min_edge_weight        : float
      selection_mode         : 'best_per_target'|'best_per_srna'|'all'
      top_n_per_group        : int  (with best_per_*, keep the top N partners per
                                     group ranked by edge weight; 0 = best only/1)
    """
    # ---- 0. Base table of sRNA-Target pairs --------------------------------
    # Prediction CSV is preferred (it carries the program scores). If it is left
    # empty, the pairs are built from the WGCNA edges file instead, and the
    # energy/probability step is skipped (no scores available).
    pred_path = str(p.get("predictions_csv", "")).strip()
    edges_path = str(p.get("module_edges_file", "")).strip()

    # ---- Optional per-step intermediate tables -----------------------------
    # When 'save_intermediate' is set, a separate CSV is written after every
    # filtering step so each step's result can be inspected on its own. The
    # files are named after the final Output CSV with a step suffix, e.g.
    # filtered_weight_step1_base.csv, _step2_energy.csv, ... in the same folder.
    save_intermediate = bool(p.get("save_intermediate", False))
    _inter_step = {"n": 0}

    def _save_step(df, label):
        """Write the current table for one filtering step (if enabled)."""
        if not save_intermediate:
            return
        _inter_step["n"] += 1
        out_csv = str(p.get("output_csv", "")).strip() or "filtered.csv"
        base, ext = os.path.splitext(out_csv)
        ext = ext or ".csv"
        step_path = f"{base}_step{_inter_step['n']}_{label}{ext}"
        try:
            df.to_csv(step_path, index=False)
            log(f"  [intermediate] step {_inter_step['n']} ({label}): "
                f"{len(df)} row(s) -> {step_path}")
        except Exception as exc:  # never let table dumping break the run
            log(f"  [intermediate] could not write {label} table: {exc}")

    if pred_path:
        _check_exists(pred_path)
        pred = pd.read_csv(pred_path)
        if not {"sRNA", "Target"}.issubset(pred.columns):
            raise ValueError("Prediction file must contain at least 'sRNA' and 'Target' columns.")
        have_scores = True
        pair_source = "Prediction CSV"
    else:
        if not edges_path:
            raise ValueError(
                "No Prediction CSV and no WGCNA edges file were given. Provide a "
                "Prediction CSV (scored sRNA-target pairs) or an edges file to "
                "build the pairs from the co-expression network."
            )
        _check_exists(edges_path)
        srna_prefix = str(p.get("srna_prefix", "sRNA")).strip()
        # Stream the (possibly huge) edges file in chunks and keep only the
        # sRNA-target edges, so we never hold the whole file in memory at once.
        parts = []
        checked_cols = False
        for chunk in pd.read_csv(edges_path, sep="\t", chunksize=250_000):
            if not checked_cols:
                if not {"fromNode", "toNode"}.issubset(chunk.columns):
                    raise ValueError("edges file must contain fromNode and toNode columns.")
                checked_cols = True
            from_s = chunk["fromNode"].astype(str).str.startswith(srna_prefix)
            to_s = chunk["toNode"].astype(str).str.startswith(srna_prefix)
            keep = from_s ^ to_s  # exactly one endpoint is an sRNA
            if not keep.any():
                continue
            sub = chunk.loc[keep]
            parts.append(pd.DataFrame({
                "sRNA": np.where(from_s.loc[keep], sub["fromNode"], sub["toNode"]),
                "Target": np.where(from_s.loc[keep], sub["toNode"], sub["fromNode"]),
            }))
        pred = (pd.concat(parts, ignore_index=True).drop_duplicates()
                if parts else pd.DataFrame(columns=["sRNA", "Target"]))
        if pred.empty:
            raise ValueError(
                f"No sRNA-target edges found: no edges have exactly one node whose "
                f"name starts with the sRNA prefix '{srna_prefix}'. Adjust the prefix."
            )
        have_scores = False
        pair_source = f"WGCNA edges file (sRNA prefix '{srna_prefix}')"

    # ---- Parameter summary (printed first, before the funnel) --------------
    deg_files = p.get("deg_files") or []
    deg_enabled = bool(p.get("deg_enabled", True)) and bool(deg_files)
    module_enabled = bool(p.get("module_filter_enabled", False))
    minw_enabled = bool(p.get("min_weight_enabled", False))
    sel_enabled = bool(p.get("selection_enabled", False))
    loc_enabled = bool(p.get("location_filter_enabled", False))
    energy_enabled = p.get("energy_enabled", {m: True for m in _ENERGY_METRICS})

    sep = "=" * 64
    log(sep)
    log("SELECTED FILTER PARAMETERS")
    log(f"  Pair source          : {pair_source}  ({len(pred)} pairs)")
    if have_scores:
        on_metrics = [m for m in _ENERGY_METRICS if energy_enabled.get(m, True)]
        emp = int(p.get("energy_min_pass", 0) or 0)
        log(f"  Energy metrics ON    : {', '.join(on_metrics) if on_metrics else 'none'}")
        log(f"  Energy consensus     : {'ALL enabled (AND)' if emp == 0 else 'at least ' + str(emp)}")
    else:
        log("  Energy metrics       : n/a (pairs built from edges)")
    if deg_enabled:
        log(f"  DEG filter           : ON  ({len(deg_files)} file(s))")
        log(f"    padj cutoff        : <= {p.get('padj_cutoff', 0.05)}")
        log(f"    baseMean           : {'>= ' + str(p.get('deg_min_basemean')) if p.get('deg_basemean_enabled') else 'OFF'}")
        band = []
        if p.get("lfc_up_enabled"):
            band.append(f">= {p.get('lfc_up')}")
        if p.get("lfc_down_enabled"):
            band.append(f"<= {p.get('lfc_down')}")
        log(f"    log2FC band        : {' or '.join(band) if band else 'OFF (keep any log2FC)'}")
        log(f"    consensus type     : {p.get('consensus_mode', 'sign')} in >= {p.get('min_strains_consistent', 1)} file(s)")
        dirs = []
        if p.get("srna_dir_enabled") and p.get("srna_direction", "any") != "any":
            dirs.append(f"sRNA={p.get('srna_direction')}")
        if p.get("target_dir_enabled") and p.get("target_direction", "any") != "any":
            dirs.append(f"target={p.get('target_direction')}")
        if p.get("pair_rel_enabled") and p.get("pair_relationship", "any") != "any":
            dirs.append(f"pair={p.get('pair_relationship')}")
        log(f"    direction filters  : {', '.join(dirs) if dirs else 'OFF'}")
    else:
        log("  DEG filter           : OFF")
    log(f"  Module filter        : {p.get('module_mode') if module_enabled else 'OFF'}")
    log(f"  Require network edge : {bool(p.get('require_in_network', True))}")
    log(f"  Min edge weight      : {'>= ' + str(p.get('min_edge_weight')) if minw_enabled else 'OFF'}")
    if sel_enabled:
        tn = int(p.get("top_n_per_group", 0) or 0)
        log(f"  Best-pair selection  : {p.get('selection_mode')} "
            f"({'top ' + str(tn) if tn > 0 else 'single best'} per group)")
    else:
        log("  Best-pair selection  : OFF (keep all)")
    log(f"  Position filter (LAST): {', '.join(p.get('locations_keep', [])) if loc_enabled else 'OFF'}")
    log(sep)
    log(f"Total pairs: {len(pred)}")
    _save_step(pred, "base")

    # ---- 1. Energy / probability filter ------------------------------------
    log("-" * 64)
    log("[1] ENERGY / PROBABILITY")
    if have_scores:
        thr = p.get("filters_outlier", {})
        pass_cols = []
        used = []
        for metric, (direction, default_thr) in _ENERGY_METRICS.items():
            if not energy_enabled.get(metric, True):
                continue
            if metric not in pred.columns:
                log(f"  [skip] {metric}: not present in the Prediction file.")
                continue
            t = thr.get(metric + ("_max" if direction == "max" else "_min"),
                        thr.get(metric, default_thr))
            col = _numeric_col(pred[metric])
            cond = (col <= t) if direction == "max" else (col >= t)
            pass_cols.append(cond)
            used.append(metric)
            sym = "<=" if direction == "max" else ">="
            log(f"  {metric} {sym} {t}: {int((~cond).sum())} pair(s) fail this metric")
        if pass_cols:
            passes = pd.concat(pass_cols, axis=1).sum(axis=1)
            n_enabled = len(pass_cols)
            min_pass = int(p.get("energy_min_pass", 0) or 0) or n_enabled
            min_pass = max(1, min(min_pass, n_enabled))
            before = len(pred)
            pred = pred.loc[passes >= min_pass].copy()
            log(f"  => consensus >= {min_pass}/{n_enabled} metrics: "
                f"removed {before - len(pred)}, kept {len(pred)}")
        else:
            log("  Energy filter skipped (no enabled metric columns present).")
    else:
        log("  Skipped (pairs built from edges; no scores available).")
    _save_step(pred, "energy")

    # ---- 2. DEG (padj -> baseMean -> log2FC -> consensus) ------------------
    log("-" * 64)
    log("[2] DIFFERENTIAL EXPRESSION (DEG)")
    if deg_enabled:
        info = _build_deg_status(
            deg_files,
            p.get("padj_cutoff", 0.05),
            p.get("min_strains_consistent", 1),
            log,
            min_basemean=p.get("deg_min_basemean", 0.0),
            basemean_enabled=bool(p.get("deg_basemean_enabled", False)),
            lfc_up_enabled=bool(p.get("lfc_up_enabled", False)),
            lfc_up=float(p.get("lfc_up", 1.0)),
            lfc_down_enabled=bool(p.get("lfc_down_enabled", False)),
            lfc_down=float(p.get("lfc_down", -1.0)),
            consensus_mode=p.get("consensus_mode", "sign"),
        )
        status_map = info["status_map"]
        lfc_map = info["log2fc_map"]
        pred["sRNA_DEG"] = pred["sRNA"].map(status_map)
        pred["Target_DEG"] = pred["Target"].map(status_map)
        pred["sRNA_log2FC"] = pred["sRNA"].map(lfc_map)
        pred["Target_log2FC"] = pred["Target"].map(lfc_map)

        # Pair funnel: BOTH endpoints must survive each cascade level. Each gene
        # set is nested in the previous, so the pair counts fall monotonically.
        cur = len(pred)
        m = _both_in(pred, info["genes_padj"]); after = int(m.sum())
        log(f"  padj (<= {p.get('padj_cutoff', 0.05)}): removed {cur - after} -> {after}")
        pred = pred[m]; cur = after
        m = _both_in(pred, info["genes_basemean"]); after = int(m.sum())
        log(f"  baseMean: removed {cur - after} -> {after}")
        pred = pred[m]; cur = after
        m = _both_in(pred, info["genes_log2fc"]); after = int(m.sum())
        log(f"  log2FC band: removed {cur - after} -> {after}")
        pred = pred[m]; cur = after
        m = _both_in(pred, set(status_map)); after = int(m.sum())
        log(f"  consensus [{p.get('consensus_mode', 'sign')} in >= "
            f"{p.get('min_strains_consistent', 1)} file(s)]: removed {cur - after} -> {after}")
        pred = pred[m].copy(); cur = after

        # Optional direction / relationship sub-filters (after consensus).
        if p.get("srna_dir_enabled") and p.get("srna_direction", "any") != "any":
            sdir = p["srna_direction"]
            before = len(pred); pred = pred[pred["sRNA_DEG"] == sdir]
            log(f"  sRNA direction = {sdir}: removed {before - len(pred)} -> {len(pred)}")
        if p.get("target_dir_enabled") and p.get("target_direction", "any") != "any":
            tdir = p["target_direction"]
            before = len(pred); pred = pred[pred["Target_DEG"] == tdir]
            log(f"  target direction = {tdir}: removed {before - len(pred)} -> {len(pred)}")
        if p.get("pair_rel_enabled") and p.get("pair_relationship", "any") != "any":
            rel = p["pair_relationship"]
            before = len(pred)
            if rel == "same":
                pred = pred[pred["sRNA_DEG"] == pred["Target_DEG"]]
            elif rel == "opposite":
                pred = pred[pred["sRNA_DEG"].notna() & pred["Target_DEG"].notna()
                            & (pred["sRNA_DEG"] != pred["Target_DEG"])]
            log(f"  sRNA vs target = {rel}: removed {before - len(pred)} -> {len(pred)}")
        log(f"  => DEG kept: {len(pred)}")
    else:
        for c in ("sRNA_DEG", "Target_DEG", "sRNA_log2FC", "Target_log2FC"):
            pred[c] = np.nan
        log("  Skipped (DEG filter off or no DEG files).")
    _save_step(pred, "deg")

    # ---- 3. WGCNA (module -> network edge -> min weight -> best pair) ------
    log("-" * 64)
    log("[3] WGCNA (module & network)")

    # Merge the edge weight onto every pair. The edges file can have millions of
    # rows, so we STREAM it in chunks and keep only the weights for the pairs we
    # actually still have, instead of loading and keying the whole file at once.
    # Loading it whole can exhaust memory and make the OS kill the app during the
    # WGCNA stage (a hard crash Python cannot catch).
    pred = pred.copy()
    if edges_path:
        _check_exists(edges_path)
        pred["pair"] = _pair_key(pred["sRNA"], pred["Target"]).to_numpy()
        wanted = set(pred["pair"].tolist())
        best_weight = {}
        n_rows = 0
        checked_cols = False
        for chunk in pd.read_csv(edges_path, sep="\t", chunksize=250_000):
            if not checked_cols:
                if not {"fromNode", "toNode", "weight"}.issubset(chunk.columns):
                    raise ValueError("edges file must contain fromNode, toNode, weight")
                checked_cols = True
            n_rows += len(chunk)
            keys = _pair_key(chunk["fromNode"], chunk["toNode"]).to_numpy()
            wts = _numeric_col(chunk["weight"]).to_numpy()
            for k, wv in zip(keys, wts):
                if k in wanted and wv == wv:  # wv == wv skips NaN
                    prev = best_weight.get(k)
                    if prev is None or wv > prev:  # keep the strongest edge
                        best_weight[k] = wv
        pred["weight"] = pred["pair"].map(best_weight).astype(float)
        pred = pred.drop(columns=["pair"])
        log(f"  scanned {n_rows} edge(s); matched {len(best_weight)} of "
            f"{len(wanted)} pair(s) to a network weight.")
    else:
        pred["weight"] = np.nan

    # 3a. Module color relationship
    if module_enabled and p.get("module_nodes_file"):
        _check_exists(p["module_nodes_file"])
        nodes = pd.read_csv(p["module_nodes_file"], sep="\t")
        nodes = nodes.rename(columns={"nodeAttr[nodesPresent, ]": "nodeAttr"})
        if not {"nodeName", "nodeAttr"}.issubset(nodes.columns):
            raise ValueError("nodes file must contain nodeName, nodeAttr")
        nm = nodes.set_index("nodeName")["nodeAttr"].to_dict()
        pred["sRNA_Module"] = pred["sRNA"].map(nm)
        pred["Target_Module"] = pred["Target"].map(nm)
        both = pred["sRNA_Module"].notna() & pred["Target_Module"].notna()
        module_mode = p.get("module_mode", "same")
        before = len(pred)
        if module_mode == "same":
            pred = pred[both & (pred["sRNA_Module"] == pred["Target_Module"])]
        elif module_mode == "different":
            pred = pred[both & (pred["sRNA_Module"] != pred["Target_Module"])]
        log(f"  module ({module_mode} color): removed {before - len(pred)} -> {len(pred)}")
    elif module_enabled:
        log("  module filter enabled but no nodes file given -> skipped.")
    else:
        log("  module filter: OFF")

    # 3b. Require a network edge
    if p.get("require_in_network", True):
        before = len(pred)
        pred = pred[pred["weight"].notna()]
        log(f"  require network edge: removed {before - len(pred)} -> {len(pred)}")
    else:
        log("  require network edge: OFF")

    # 3c. Minimum edge weight
    if minw_enabled:
        min_w = float(p.get("min_edge_weight", 0.0) or 0.0)
        before = len(pred)
        pred = pred[(pred["weight"].isna()) | (pred["weight"] >= min_w)]
        log(f"  min edge weight (>= {min_w}): removed {before - len(pred)} -> {len(pred)}")
    else:
        log("  min edge weight: OFF")
    # Save the result of the module + network-edge requirements on its own,
    # BEFORE the best-pair-per-target prioritization, so this table reflects
    # every pair that passed same-module / require-edge / min-weight (it is not
    # yet reduced to one pair per target).
    _save_step(pred, "wgcna")

    # 3d. Keep best pair per group (best_per_srna / best_per_target) + Top N
    if sel_enabled:
        selection = p.get("selection_mode", "best_per_target")
        key = "Target" if selection == "best_per_target" else "sRNA"
        top_n = int(p.get("top_n_per_group", 0) or 0)
        before = len(pred)
        pred = pred.sort_values([key, "weight"], ascending=[True, False])
        if top_n > 0:
            pred = pred.groupby(key, sort=False, group_keys=False).head(top_n)
            log(f"  {selection} top {top_n} per {key}: "
                f"removed {before - len(pred)} -> {len(pred)}")
        else:
            pred = pred.drop_duplicates(key)
            log(f"  {selection} single best per {key}: "
                f"removed {before - len(pred)} -> {len(pred)}")
    else:
        log("  best-pair selection: OFF (keep all)")
    # Save the result AFTER best-pair-per-target prioritization on its own
    # table, so the reduction caused only by keeping the single best partner
    # per group is visible separately from the module/network requirements.
    _save_step(pred, "bestpair")

    # ---- 4. Position classification (ALWAYS LAST) --------------------------
    log("-" * 64)
    log("[4] POSITION CLASSIFICATION (last step)")
    if loc_enabled:
        ann_path = str(p.get("srna_annotation_file", "")).strip()
        keep = p.get("locations_keep", [])
        if not ann_path:
            raise ValueError(
                "Location filter is enabled but no annotation file was given. "
                "Provide the sRNA annotation table from the 01 Position "
                "Classification tab, or disable the location filter."
            )
        if not keep:
            raise ValueError(
                "Location filter is enabled but no location categories are "
                "selected. Tick at least one location to keep."
            )
        _check_exists(ann_path)
        if ann_path.lower().endswith((".xlsx", ".xls")):
            ann = pd.read_excel(ann_path)
        elif ann_path.lower().endswith(".tsv"):
            ann = pd.read_csv(ann_path, sep="\t")
        else:
            ann = pd.read_csv(ann_path)
        id_col = str(p.get("location_id_col", "exon_id")).strip() or "exon_id"
        if id_col not in ann.columns or "location" not in ann.columns:
            raise ValueError(
                f"Annotation file must contain '{id_col}' and 'location' columns "
                f"(found: {list(ann.columns)})."
            )
        loc_map = ann.set_index(ann[id_col].astype(str))["location"].to_dict()
        pred["sRNA_Location"] = pred["sRNA"].astype(str).map(loc_map)
        counts = pred["sRNA_Location"].value_counts(dropna=False)
        log("  current pairs by sRNA position:")
        for cat, cnt in counts.items():
            label = "(not in annotation)" if pd.isna(cat) else cat
            log(f"    {label}: {cnt}")
        keep_set = set(keep)
        before = len(pred)
        pred = pred[pred["sRNA_Location"].isin(keep_set)].copy()
        log(f"  keep {', '.join(keep)}: removed {before - len(pred)} -> {len(pred)}")
    else:
        log("  Skipped.")
    _save_step(pred, "position")

    log(sep)
    out = p["output_csv"]
    pred.to_csv(out, index=False)
    log(f"FINAL interactions: {len(pred)}")
    log(f"Saved: {out}")
    return out


# =====================================================================
# Combine prediction-program outputs into a single Prediction table
# =====================================================================
#
# Reads the native output of up to four sRNA-target prediction programs and
# merges them, on the (Target, sRNA) pair, into one table with the same columns
# as Prediction_test.csv:
#
#   Target, sRNA, E_intaRNA, p_intaRNA, E_Rnaplex, p_Rnaplex,
#   E_TargetRNA3, p_TargetRNA3, Probability_TargetRNA3, Probability_sRNARFTarget
#
# intaRNA and RNAplex only report an interaction energy, so their p-value column
# is computed as an empirical p-value (the ECDF of the energies: the fraction of
# pairs whose energy is at least as low). Any program left blank is skipped.

_COMBINE_OUTPUT_ORDER = [
    "Target", "sRNA",
    "E_intaRNA", "p_intaRNA",
    "E_Rnaplex", "p_Rnaplex",
    "E_TargetRNA3", "p_TargetRNA3", "Probability_TargetRNA3",
    "Probability_sRNARFTarget",
]


def _find_col(cols, *candidates):
    low = {str(c).lower().strip(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def parse_intarna(path, log):
    """intaRNA ';'-separated table: id1=Target, id2=sRNA, E=interaction energy."""
    _check_exists(path)
    df = pd.read_csv(path, sep=";")
    tcol = _find_col(df.columns, "id1", "Target", "mRNA", "mRNA_ID")
    scol = _find_col(df.columns, "id2", "sRNA", "sRNA_ID")
    ecol = _find_col(df.columns, "E", "E_intaRNA", "Energy")
    if not (tcol and scol and ecol):
        raise ValueError(
            f"intaRNA file must have target/sRNA/energy columns; found {list(df.columns)}")
    out = pd.DataFrame({
        "Target": df[tcol].astype(str),
        "sRNA": df[scol].astype(str),
        "E_intaRNA": pd.to_numeric(df[ecol], errors="coerce"),
    })
    pcol = _find_col(df.columns, "p_intaRNA", "p-value", "pvalue", "p")
    if pcol:
        out["p_intaRNA"] = pd.to_numeric(df[pcol], errors="coerce")
    out = out.sort_values("E_intaRNA").drop_duplicates(["Target", "sRNA"])
    log(f"  intaRNA: {len(out)} unique target-sRNA pairs")
    return out


def parse_rnaplex(path, log):
    """RNAplex output: blocks of '>target', '>srna', structure line w/ '(energy)'."""
    _check_exists(path)
    energy_re = re.compile(r"\(\s*(-?\d+\.?\d*)\s*\)")
    rows = []
    target = srna = None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if target is None:
                    target = line[1:].strip()
                else:
                    srna = line[1:].strip()
            elif line and target is not None and srna is not None:
                m = energy_re.search(line)
                if m:
                    rows.append((target, srna, float(m.group(1))))
                target = srna = None
    out = pd.DataFrame(rows, columns=["Target", "sRNA", "E_Rnaplex"])
    out = out.sort_values("E_Rnaplex").drop_duplicates(["Target", "sRNA"])
    log(f"  RNAplex: {len(out)} unique target-sRNA pairs")
    return out


def parse_targetrna3(path, log):
    """TargetRNA3 tab table. Needs an sRNA column (plus Target) to be combinable."""
    _check_exists(path)
    df = pd.read_csv(path, sep="\t")
    scol = _find_col(df.columns, "sRNA", "sRNA_ID", "sRNA_id")
    tcol = _find_col(df.columns, "Target", "mRNA", "mRNA_ID")
    ecol = _find_col(df.columns, "Energy", "E_TargetRNA3", "E")
    pcol = _find_col(df.columns, "P-value", "p_TargetRNA3", "pvalue", "p")
    probcol = _find_col(df.columns, "Probability", "Probability_TargetRNA3")
    if not (scol and tcol):
        raise ValueError(
            "TargetRNA3 file must contain an 'sRNA' column and a 'Target' column to "
            f"be combinable; found {list(df.columns)}. (TargetRNA3 runs one sRNA at "
            "a time; add an sRNA column identifying which sRNA each block belongs to.)")
    out = pd.DataFrame({"Target": df[tcol].astype(str), "sRNA": df[scol].astype(str)})
    if ecol:
        out["E_TargetRNA3"] = pd.to_numeric(df[ecol], errors="coerce")
    if pcol:
        out["p_TargetRNA3"] = pd.to_numeric(df[pcol], errors="coerce")
    if probcol:
        out["Probability_TargetRNA3"] = pd.to_numeric(df[probcol], errors="coerce")
    if ecol:
        out = out.sort_values("E_TargetRNA3")
    out = out.drop_duplicates(["Target", "sRNA"])
    log(f"  TargetRNA3: {len(out)} unique target-sRNA pairs")
    return out


def parse_srnarftarget(path, log):
    """sRNARFTarget table (tab or comma): sRNA_ID, mRNA_ID, Prediction_Probability."""
    _check_exists(path)
    df = pd.read_csv(path, sep="\t")
    if df.shape[1] == 1:  # fall back to comma despite a .csv extension
        df = pd.read_csv(path)
    scol = _find_col(df.columns, "sRNA_ID", "sRNA")
    tcol = _find_col(df.columns, "mRNA_ID", "Target", "mRNA")
    pcol = _find_col(df.columns, "Prediction_Probability", "Probability_sRNARFTarget", "Probability")
    if not (scol and tcol and pcol):
        raise ValueError(
            f"sRNARFTarget file needs sRNA/mRNA/Probability columns; found {list(df.columns)}")
    out = pd.DataFrame({
        "Target": df[tcol].astype(str),
        "sRNA": df[scol].astype(str),
        "Probability_sRNARFTarget": pd.to_numeric(df[pcol], errors="coerce"),
    })
    out = out.sort_values("Probability_sRNARFTarget", ascending=False).drop_duplicates(["Target", "sRNA"])
    log(f"  sRNARFTarget: {len(out)} unique target-sRNA pairs")
    return out


def _empirical_p(energy):
    """Empirical p-value = ECDF of the energies (fraction of pairs <= this one)."""
    s = pd.to_numeric(energy, errors="coerce")
    n = s.notna().sum()
    return s.rank(method="average") / n if n else s


def run_combine_predictions(p, log):
    """Merge up to four prediction-program outputs into a Prediction_test-style table.

    Recognized params (any file path may be blank to skip that program):
      intarna_file, rnaplex_file, targetrna3_file, srnarftarget_file
      compute_empirical_p   : bool  (fill p_intaRNA / p_Rnaplex from the energy ECDF)
      require_all_programs  : bool  (inner join: keep only pairs predicted by every
                                     provided program; default False = outer join)
      relabel_ids           : bool  (rename Target->Gene#, sRNA->sRNA#)
      mapping_csv           : str   (optional path to save the relabel mapping)
      output_csv            : str   (required output path)
    """
    parsers = [
        ("intarna_file", parse_intarna),
        ("rnaplex_file", parse_rnaplex),
        ("targetrna3_file", parse_targetrna3),
        ("srnarftarget_file", parse_srnarftarget),
    ]
    tables = []
    for key, fn in parsers:
        path = str(p.get(key, "")).strip()
        if path:
            tables.append(fn(path, log))
    if not tables:
        raise ValueError(
            "No prediction files were given. Provide at least one program output "
            "(intaRNA, RNAplex, TargetRNA3 or sRNARFTarget).")

    how = "inner" if p.get("require_all_programs", False) and len(tables) > 1 else "outer"
    merged = tables[0]
    for t in tables[1:]:
        merged = pd.merge(merged, t, on=["Target", "sRNA"], how=how)
    log(f"Merged ({how} join): {len(merged)} target-sRNA pairs")

    if p.get("compute_empirical_p", True):
        if "E_intaRNA" in merged.columns and "p_intaRNA" not in merged.columns:
            merged["p_intaRNA"] = _empirical_p(merged["E_intaRNA"])
        if "E_Rnaplex" in merged.columns and "p_Rnaplex" not in merged.columns:
            merged["p_Rnaplex"] = _empirical_p(merged["E_Rnaplex"])

    if p.get("relabel_ids", False):
        tmap = {t: f"Gene{i+1}" for i, t in enumerate(sorted(merged["Target"].unique()))}
        smap = {s: f"sRNA{i+1}" for i, s in enumerate(sorted(merged["sRNA"].unique()))}
        mp = str(p.get("mapping_csv", "")).strip()
        if mp:
            pd.DataFrame({
                "type": ["Target"] * len(tmap) + ["sRNA"] * len(smap),
                "original": list(tmap) + list(smap),
                "renamed": list(tmap.values()) + list(smap.values()),
            }).to_csv(mp, index=False)
            log(f"Saved ID mapping: {mp}")
        merged["Target"] = merged["Target"].map(tmap)
        merged["sRNA"] = merged["sRNA"].map(smap)

    cols = [c for c in _COMBINE_OUTPUT_ORDER if c in merged.columns]
    merged = merged[cols]
    out = p["output_csv"]
    out_dir = os.path.dirname(out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    merged.to_csv(out, index=False)
    log(f"Combined prediction table: {len(merged)} rows, {len(cols)} columns.")
    log(f"Saved: {out}")
    return out
