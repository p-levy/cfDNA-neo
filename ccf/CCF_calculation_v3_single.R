#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  suppressWarnings({
    library(data.table)
    library(tidyverse)
    library(bedtoolsr)
	library(argparse)
  })
})

# Function to standardize chromosome names
standardize_chr <- function(chr_vector) {
  # Convert to character first
  chr_vector <- as.character(chr_vector)
  
  # Remove "chr" prefix if present
  chr_vector <- gsub("^chr", "", chr_vector, ignore.case = TRUE)
  
  return(chr_vector)
}

# Argument parser
parser <- ArgumentParser(description = "Estimate cancer cell fraction (CCF) from variant allele frequencies, purity and CNA profiles")

parser$add_argument("--patient", type = "character", required = TRUE,
                    help = "Patient ID") #  patient <- "test"
parser$add_argument("--variants_counts_path", type = "character", required = TRUE,
                    help = "Path to variant read counts TSV") # variants_counts_path <- "../test/variant-counts/test_variants_counts.tsv"
parser$add_argument("--outdir", type = "character", required = TRUE,
                    help = "Output directory") # outdir <- "../test/ccf"
parser$add_argument("--segs", type = "character", required = TRUE,
                    help = "ASCAT raw segments") # segs <- fread("../test/cna/test_FrTu.hg38.segments_raw.txt")
parser$add_argument("--purity", type = "character", required = TRUE,
                    help = "Purity value for Tumor sample (ASCAT or 'vaf')") # purity <- 0.41 # Use ASCAT purity or Write "vaf" if peak_vaf has to be used!

# Filters
parser$add_argument("--min_tvaf", type = "double", default = 0.03,
                    help = "Minimum tumor VAF (default: 0.03)") # min_tvaf <- 0.03
parser$add_argument("--min_alt", type = "integer", default = 4,
                    help = "Minimum ALT reads (default: 4)") # min_alt <- 4
parser$add_argument("--min_cov", type = "integer", default = 9,
                    help = "Minimum coverage (default: 9)") # min_cov <- 9

# Optional inputs
parser$add_argument("--bed_exome", type = "character", required = FALSE,
                    help = "Exome capture BED file") # bed_exome <- fread("../test/bed/test.hg38.bed")
parser$add_argument("--nsm_annot", type = "character", required = FALSE,
                    help = "NSM annotation table")  # nsm_annot <- fread("../test/annot/test_nsm_annot.tsv")
parser$add_argument("--neo", type = "character", required = FALSE,
                    help = "Neo output (WiGiTS)")  # neo <- fread("../test/neo/test_neo.neoepitope.tsv")

# Load args
args <- parser$parse_args()

patient <- args$patient
variants_counts_path <- args$variants_counts_path
outdir <- args$outdir
segs  <- fread(args$segs)
purity <- args$purity

if (!is.null(args$bed_exome)) {bed_exome <- fread(args$bed_exome)}
if (!is.null(args$nsm_annot) & !is.null(args$neo)) {stop("Only provide one NSM annotation file (--nsm_annot OR --neo)")}
if (!is.null(args$nsm_annot)) {nsm_annot <- fread(args$nsm_annot)}
if (!is.null(args$neo)) {neo <- fread(args$neo)}

# Filters
min_tvaf <- args$min_tvaf
min_alt <- args$min_alt
min_cov <- args$min_cov
min_callers_snv <- args$min_callers_snv
min_callers_indels <- args$min_callers_indels

# Create output dir if not existing already
system(paste0("mkdir -p ", outdir))

# Load mutation count data
variants_counts <- fread(variants_counts_path)

# Standardize chromosome names in variants_counts
variants_counts$CHROM <- standardize_chr(variants_counts$CHROM)

# Standardize chromosome names in segs (assuming chromosome is in column 2)
segs$chr <- standardize_chr(segs$chr)

# OPTIONNAL: Add NSM info (+ wt/mut epitope sequences, immunogenicity, etc.)
if (!is.null(args$nsm_annot)) {
    # Standardize chromosome names in nsm_annot
    nsm_annot$CHROM <- standardize_chr(nsm_annot$CHROM)
    variants_counts <- variants_counts %>% left_join(nsm_annot %>% mutate(coding_consequence = "non_synonymous"), by = c("CHROM", "POS", "REF", "ALT"))
    variants_counts$coding_consequence <- replace_na(variants_counts$coding_consequence, "synonymous_or_noncoding")
}

# OPTIONNAL: Add info from Neo output if provided
if (!is.null(args$neo)) {
	variants_counts <- neo %>% full_join(variants_counts %>% mutate(coding_consequence = "non_synonymous"), by = "VariantInfo")
	variants_counts$coding_consequence <- replace_na(variants_counts$coding_consequence, "synonymous_or_noncoding")
}

# Apply filters (on tvaf, min_cov, min_callers...)
variants_counts <- variants_counts %>%
    filter(vaf_Tumor >= min_tvaf) %>%
    filter(ALT_counts_Tumor >= min_alt) %>%
    filter(cov_Tumor >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | vaf_Tumor >= 5 * vaf_Normal)

# OPTIONNAL: Use bedtoolsr to filter the mutations to keep only the ones falling in the exome bed file
if (!is.null(args$bed_exome)) {
    # Standardize chromosome names in bed_exome
    bed_exome[[1]] <- standardize_chr(bed_exome[[1]])
    
    bed_mut <- variants_counts %>%
        dplyr::select(CHROM, POS) %>%
        mutate(end = POS) %>%
        dplyr::rename(start = POS, chrom = CHROM)
    bed_mut_exome <- bt.intersect(bed_mut, bed_exome, u = TRUE) # u=TRUE to only report one entry per mutation, if at least it is in one of the bed intervals
    bed_mut_exome <- bed_mut_exome %>% dplyr::rename(CHROM = V1, POS = V2)
    bed_mut_exome$CHROM <- standardize_chr(bed_mut_exome$CHROM)
    variants_counts <- bed_mut_exome %>%
        dplyr::select(CHROM, POS) %>%
        left_join(variants_counts, by = c("CHROM", "POS"))

	# Get the copy number from the segment files and annotate the mutation table
	bed_seg_Tumor <- segs %>% dplyr::select(2:6)
	intersect_Tumor <- bt.intersect(bed_mut_exome, bed_seg_Tumor, wb = T)
	intersect_Tumor <- intersect_Tumor %>%
		dplyr::select(1, 2, 7, 8) %>%
		dplyr::rename(CHROM = V1, POS = V2, nMajor = V7, nMinor = V8)
} else {
	# Get the copy number from the segment files and annotate the mutation table
	bed_mut <- variants_counts %>%
        dplyr::select(CHROM, POS) %>%
        mutate(end = POS) %>%
        dplyr::rename(start = POS, chrom = CHROM)
	bed_seg_Tumor <- segs %>% dplyr::select(2:6)
	intersect_Tumor <- bt.intersect(bed_mut, bed_seg_Tumor, wb = T)
	intersect_Tumor <- intersect_Tumor %>%
		dplyr::select(1, 2, 7, 8) %>%
		dplyr::rename(CHROM = V1, POS = V2, nMajor = V7, nMinor = V8)
}

# Final mutation table with CN annotated
intersect_Tumor$CHROM <- standardize_chr(intersect_Tumor$CHROM)
variants_counts <- variants_counts %>%
    left_join(intersect_Tumor, by = c("CHROM", "POS"))

# Add CCF calculation
if (str_to_lower(purity) == "vaf") {
    # Calculate purity based on VAF of most abundant cluster
    table_Tumor <- variants_counts
    # Find x value for peak density
    peak_vaf_Tumor <- density(table_Tumor$vaf_Tumor)$x[which.max(density(table_Tumor$vaf_Tumor)$y)]
    purity <- 2 * peak_vaf_Tumor
} else {
    purity <- as.numeric(purity)
}

# CCF = (VAF / Purity)*((purity * Tumour_TotalCN) + (Normal_TotalCN*( 1 – Purity ))))
ccf <- variants_counts %>%
    mutate(ccf = (vaf_Tumor / purity) * ((1 - purity) * 2 + purity * (nMajor  + nMinor)),
	purity = purity)

ccf <- ccf %>% distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)

# Export table
write.table(ccf, file = paste0(outdir, "/", patient, "_all_mutations_CCF.tsv"), row.names = F, quote = F, sep = "\t")
