
import pandas as pd
import numpy as np
import re

###############################################
# CONFIGURATION
###############################################

GTF_FILE = "BMB.gtf"

GTF_FEATURE = "exon"          # feature type in GTF
GTF_ID_FIELD = "exon_id"      # attribute key to extract

SRNA_PREFIX = "srn_"
GENE_PREFIX = "SABB"

SRNA_RENAME_PREFIX = "sRNA"
GENE_RENAME_PREFIX = "Gene"

UTR_DISTANCE = 150


###############################################
# FUNCTIONS
###############################################

def read_gtf(gtf_path):

    gtf_columns = [
        "seqname","source","feature","start","end",
        "score","strand","frame","attribute"
    ]

    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=gtf_columns
    )

    return gtf


def extract_attribute(attribute, key):

    if pd.isna(attribute):
        return None

    pattern = rf'{re.escape(key)} "([^"]+)"'
    match = re.search(pattern, attribute)

    return match.group(1) if match else None


###############################################
# LOAD GTF
###############################################

gtf = read_gtf(GTF_FILE)

gtf = gtf[gtf["feature"] == GTF_FEATURE].copy()

gtf["raw_id"] = gtf["attribute"].apply(lambda x: extract_attribute(x, GTF_ID_FIELD))

gtf = gtf[["raw_id","start","end","strand"]].copy()

gtf["length"] = gtf["end"] - gtf["start"]

gtf = gtf[gtf["raw_id"].notna()].copy()


###############################################
# FILTER SRNA AND GENES
###############################################

filtered = gtf[
    gtf["raw_id"].astype(str).str.startswith((SRNA_PREFIX,GENE_PREFIX))
].copy()

srn_results = filtered[
    filtered["raw_id"].str.startswith(SRNA_PREFIX)
].copy()

gene_results = filtered[
    filtered["raw_id"].str.startswith(GENE_PREFIX)
].copy()

srn_results = srn_results.reset_index(drop=True)
gene_results = gene_results.reset_index(drop=True)


###############################################
# RENAME IDS
###############################################

srn_results["exon_id"] = [
    f"{SRNA_RENAME_PREFIX}{i+1}" for i in range(len(srn_results))
]

gene_results["gene_id"] = [
    f"{GENE_RENAME_PREFIX}{i+1}" for i in range(len(gene_results))
]

srn_results["original_id"] = srn_results["raw_id"]
gene_results["original_id"] = gene_results["raw_id"]

srn_results = srn_results[
    ["exon_id","original_id","start","end","strand","length"]
].copy()

gene_results = gene_results[
    ["gene_id","original_id","start","end","strand","length"]
].copy()


###############################################
# CLASSIFICATION
###############################################

srn_results["location"] = None
srn_results["gene_id"] = None

srn_results[["start","end"]] = srn_results[["start","end"]].astype(int)
gene_results[["start","end"]] = gene_results[["start","end"]].astype(int)


###############################################
# INTRAGENIC / ANTISENSE
###############################################

for _, gene in gene_results.iterrows():

    srn_len = srn_results["end"] - srn_results["start"]

    overlap_start = srn_results["start"].clip(lower=gene["start"])
    overlap_end = srn_results["end"].clip(upper=gene["end"])

    overlap_len = (overlap_end - overlap_start).clip(lower=0)

    overlap_ratio = overlap_len / srn_len

    mask = overlap_ratio >= 0.9

    same_strand = srn_results["strand"] == gene["strand"]

    srn_results.loc[mask & same_strand, "location"] = "Intragenic"
    srn_results.loc[mask & ~same_strand, "location"] = "Antisense"

    srn_results.loc[mask, "gene_id"] = gene["gene_id"]


###############################################
# UTR CLASSIFICATION
###############################################

unclassified = srn_results["location"].isna()

for _, gene in gene_results.iterrows():

    subset = srn_results[unclassified].copy()

    if subset.empty:
        break

    same_strand = subset["strand"] == gene["strand"]

    subset = subset[same_strand]

    srn_len = subset["end"] - subset["start"]

    upstream_start = subset["start"].clip(lower=gene["start"] - UTR_DISTANCE)
    upstream_end = subset["end"].clip(upper=gene["start"])

    upstream_len = (upstream_end - upstream_start).clip(lower=0)
    upstream_ratio = upstream_len / srn_len

    downstream_start = subset["start"].clip(lower=gene["end"])
    downstream_end = subset["end"].clip(upper=gene["end"] + UTR_DISTANCE)

    downstream_len = (downstream_end - downstream_start).clip(lower=0)
    downstream_ratio = downstream_len / srn_len


    upstream_mask = upstream_ratio >= 0.9
    downstream_mask = downstream_ratio >= 0.9


    if gene["strand"] == "+":

        srn_results.loc[subset.index[upstream_mask],"location"] = "5' UTR"
        srn_results.loc[subset.index[downstream_mask],"location"] = "3' UTR"

    else:

        srn_results.loc[subset.index[upstream_mask],"location"] = "3' UTR"
        srn_results.loc[subset.index[downstream_mask],"location"] = "5' UTR"


    srn_results.loc[
        subset.index[upstream_mask | downstream_mask],
        "gene_id"
    ] = gene["gene_id"]


###############################################
# INTERGENIC
###############################################

srn_results.loc[srn_results["location"].isna(),"location"] = "Intergenic"
srn_results.loc[srn_results["location"]=="Intergenic","gene_id"] = None


###############################################
# EXPORT
###############################################

srn_results.to_excel("sRNA_annotation.xlsx",index=False)

print("Finished.")
print(srn_results.head())
