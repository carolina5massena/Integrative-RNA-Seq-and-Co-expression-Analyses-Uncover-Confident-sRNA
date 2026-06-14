###############################################
# AUXILIARY COMMANDS – RNA-seq PIPELINE HELP
#
# This script shows:
# 1) How to download raw reads from SRA
# 2) How to convert SRA to FASTQ
# 3) How to run nf-core/rnaseq
# 4) How to run nf-core/differentialabundance
# 5) How to run the sRNA-mRNA target-prediction programs:
#       IntaRNA, RNAplex, TargetRNA3, sRNARFTarget
#
# Requirements:
# - SRA Toolkit (prefetch, fasterq-dump)
# - Nextflow
# - Docker
# - nf-core pipelines
# - IntaRNA            (https://github.com/BackofenLab/IntaRNA)
# - ViennaRNA / RNAplex (https://www.tbi.univie.ac.at/RNA/)
# - TargetRNA3         (https://github.com/btjaden/TargetRNA3)
# - sRNARFTarget       (https://github.com/BioinformaticsLabAtMUN/sRNARFTarget)
###############################################


###############################################
# STEP 1 — Download SRA file
# This retrieves SRR28417311 from NCBI SRA
###############################################

prefetch SRR28417311


###############################################
# STEP 2 — Convert SRA to FASTQ
# --split-files : separates paired-end reads
# --threads 8   : use 8 CPU threads
# -O            : output directory
###############################################

fasterq-dump SRR28417311 -O ./fastq_files --split-files --threads 8


###############################################
# nf-core RNA-seq PIPELINE
###############################################

# Repository options:
# HTTPS: https://github.com/nf-core/rnaseq.git
# SSH:   git@github.com:nf-core/rnaseq.git


###############################################
# STEP 3 — Run nf-core/rnaseq
#
# --input  : samplesheet.csv with FASTQ paths
# --outdir : output directory
# --gtf    : gene annotation file
# --fasta  : reference genome
# -profile docker : run everything inside Docker
###############################################

nextflow run nf-core/rnaseq \
    --input <samplesheet.csv> \
    --outdir <OUTDIR> \
    --gtf <genes.gtf> \
    --fasta <genome.fa> \
    -profile docker


###############################################
# nf-core DIFFERENTIAL ABUNDANCE PIPELINE
###############################################

# Repository options:
# HTTPS: https://github.com/nf-core/differentialabundance.git
# SSH:   git@github.com:nf-core/differentialabundance.git
# GitHub CLI:
# gh repo clone nf-core/differentialabundance


###############################################
# STEP 4 — Run differential abundance analysis
#
# --input     : samplesheet.csv
# --contrasts : contrasts.csv (group comparisons)
# --matrix    : counts.tsv from RNA-seq
# --gtf       : same gene annotation
# --outdir    : output directory
# -profile rnaseq,docker : use rnaseq settings + Docker
###############################################

nextflow run nf-core/differentialabundance \
     --input <samplesheet.csv> \
     --contrasts <contrasts.csv> \
     --matrix <counts.tsv> \
     --gtf <genes.gtf> \
     --outdir <OUTDIR>  \
     -profile rnaseq,docker


###############################################################################
# sRNA-mRNA TARGET PREDICTION PROGRAMS
###############################################################################
#
# These four programs each predict sRNA-mRNA interactions. Run any subset; their
# native outputs (one file per program) are then merged on the (Target, sRNA)
# pair into the single Prediction CSV used by "05 Filters"
# (columns: Target, sRNA, E_intaRNA, p_intaRNA, E_Rnaplex, p_Rnaplex,
#  E_TargetRNA3, p_TargetRNA3, Probability_TargetRNA3, Probability_sRNARFTarget).
#
# Common inputs (prepare these first):
#   sRNA.fasta   : FASTA of the sRNA sequence(s)            -> the "query"
#   mRNA.fasta   : FASTA of the candidate target/mRNA seqs  -> the "target"
#                  (e.g. CDS + upstream region, extracted from the genome/GTF)
#
# Convention used by this pipeline: target = mRNA/gene, query = sRNA.
# Example native outputs are in the "data/" folder (intaRNA_output.txt,
# Rnaplex_output.txt, TargetRNA3_output.txt, sRNARFTarget_output.csv).


###############################################
# STEP 5a — IntaRNA  (IntaRNAsTar personality)
#
# This pipeline uses the "IntaRNAsTar" personality of IntaRNA, which provides
# parameters optimized for genome-wide sRNA-target prediction (Raden et al.,
# 2020) and already emits the CSV columns the pipeline expects.
#
# -t : target sequences (mRNA / genes)
# -q : query sequences  (sRNA)
# Output is ';'-separated:  id1;id2;start1;end1;start2;end2;E
#   (id1 = target/mRNA, id2 = query/sRNA, E = interaction energy)
###############################################

IntaRNAsTar \
    -t mRNA.fasta \
    -q sRNA.fasta \
    > intaRNA_output.txt

# IntaRNAsTar is a preset personality. The two calls below are equivalent to it:
#
#   IntaRNA --personality=IntaRNAsTar -t mRNA.fasta -q sRNA.fasta > intaRNA_output.txt
#
#   IntaRNA -t mRNA.fasta -q sRNA.fasta \
#       --intLenMax=60 --intLoopMax=8 --seedNoGU --seedMinPu=0.001 \
#       --outMinPu=0.001 --outNoLP --outNoGUend --outOverlap=Q \
#       --outMode=C --outCsvCols=id1,id2,start1,end1,start2,end2,E \
#       > intaRNA_output.txt


###############################################
# STEP 5b — RNAplex (ViennaRNA Package)
#
# -q : query sequences (sRNA)
# -t : target sequences (mRNA / genes)   (-q implies -t and vice versa)
# Results are written to stdout, one interaction block per pair:
#   >target_id
#   >query_id
#   <dot-bracket> from,to : from,to (energy)
###############################################

RNAplex \
    -q sRNA.fasta \
    -t mRNA.fasta \
    > Rnaplex_output.txt


###############################################
# STEP 5c — TargetRNA3
#
# TargetRNA3 is a Python program and predicts targets for ONE sRNA at a time
# against a prepared genome directory (genomic.fna, protein.faa,
# rna_from_genomic.fna, feature_table.txt; see the TargetRNA3 User Manual).
#
# -s    : FASTA file with the sRNA sequence
# -g    : path to the genome directory
# -o    : output file (tab-separated; default is stdout)
# -prob : probability threshold (use 0 to report ALL candidate targets)
# -pval : p-value threshold      (use 1 to report ALL candidate targets)
###############################################

python TargetRNA3.py \
    -s sRNA.fa \
    -g genome_dir/ \
    -o TargetRNA3_output.txt \
    -prob 0 \
    -pval 1

# TargetRNA3 output has no "sRNA" column (it runs one sRNA at a time). When
# predicting for several sRNAs, run it once per sRNA and concatenate the results,
# adding an "sRNA" column that identifies the sRNA for each block, e.g.:
#
#   for s in sRNAs/*.fa; do
#       name=$(basename "$s" .fa)
#       python TargetRNA3.py -s "$s" -g genome_dir/ -prob 0 -pval 1 \
#           | awk -v OFS='\t' -v n="$name" 'NR==1{print "sRNA",$0; next}{print n,$0}'
#   done > TargetRNA3_output.txt


###############################################
# STEP 5d — sRNARFTarget (Nextflow pipeline)
#
# Run from inside the cloned sRNARFTarget directory, with sRNA.fasta and
# mRNA.fasta placed in that directory.
# --s : sRNA FASTA
# --m : mRNA FASTA
# Result -> sRNARFTargetResult/Prediction_probabilities.csv
#           (columns: sRNA_ID, mRNA_ID, Prediction_Probability)
###############################################

# Recommended: with the provided Docker container
nextflow run sRNARFTarget.nf \
    --s sRNA.fasta \
    --m mRNA.fasta \
    -with-docker penacastillolab/python38env

# Or, with Python installed locally:
# nextflow run sRNARFTarget.nf --s sRNA.fasta --m mRNA.fasta

cp sRNARFTargetResult/Prediction_probabilities.csv sRNARFTarget_output.csv



