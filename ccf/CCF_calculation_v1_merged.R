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
                    help = "ASCAT raw segments") # segs <- fread("../test/cna/test_FrTu.segments_raw.txt")
parser$add_argument("--purity", type = "character", required = TRUE,
                    help = "Purity value for Tumor sample (ASCAT or 'vaf')") # purity <- 0.41 # Use ASCAT purity or Write "vaf" if peak_vaf has to be used!

# Filters
parser$add_argument("--min_tvaf", type = "double", default = 0.03,
                    help = "Minimum tumor VAF (default: 0.03)") # min_tvaf <- 0.03
parser$add_argument("--min_alt", type = "integer", default = 4,
                    help = "Minimum ALT reads (default: 4)") # min_alt <- 4
parser$add_argument("--min_cov", type = "integer", default = 9,
                    help = "Minimum coverage (default: 9)") # min_cov <- 9
parser$add_argument("--min_callers_snv", type = "integer", default = 2,
                    help = "Minimum SNV callers (default: 2)") # min_callers_snv <- 2
parser$add_argument("--min_callers_indels", type = "integer", default = 2,
                    help = "Minimum INDEL callers (default: 2)") # min_callers_indels <- 1

# Optional inputs
parser$add_argument("--bed_exome", type = "character", required = FALSE,
                    help = "Exome capture BED file") # bed_exome <- fread("../test/bed/test.bed")
parser$add_argument("--nsm_annot", type = "character", required = FALSE,
                    help = "NSM annotation table")  # nsm_annot <- fread("../test/annot/test_nsm_annot.tsv")

# Load args
args <- parser$parse_args()

patient <- args$patient
variants_counts_path <- args$variants_counts_path
outdir <- args$outdir
segs  <- fread(args$segs)
purity <- args$purity

if (!is.null(args$bed_exome)) {bed_exome <- fread(args$bed_exome)}
if (!is.null(args$nsm_annot)) {nsm_annot <- fread(args$nsm_annot)}

# Filters
min_tvaf <- args$min_tvaf
min_alt <- args$min_alt
min_cov <- args$min_cov
min_callers_snv <- args$min_callers_snv
min_callers_indels <- args$min_callers_indels

# Create output dir if not existing already
system(paste0("mkdir -p ", outdir))

# Load mutation count data
variants <- fread(variants_counts_path)

# Add VAFs
variants <- variants %>%
    mutate(vaf_Tumor = as.numeric(ALT_counts_Tumor) / (as.numeric(REF_counts_Tumor) + as.numeric(ALT_counts_Tumor))) %>%
    mutate(vaf_Normal = as.numeric(ALT_counts_Normal) / (as.numeric(REF_counts_Normal) + as.numeric(ALT_counts_Normal)))

# Add Coverage
variants <- variants %>%
    mutate(cov_Tumor = as.numeric(REF_counts_Tumor) + as.numeric(ALT_counts_Tumor)) %>%
    mutate(cov_Normal = as.numeric(REF_counts_Normal) + as.numeric(ALT_counts_Normal))

# Add n_callers
variants <- variants %>%
    rowwise() %>%
    mutate(
        n_callers = ifelse(set == "Intersection", 4, length(str_split(set, "-")[[1]]))
    )

# Standardize chromosome names in variants
variants$CHROM <- standardize_chr(variants$CHROM)

# Standardize chromosome names in segs (assuming chromosome is in column 2)
segs$chr <- standardize_chr(segs$chr)

# Apply filters (on tvaf, min_cov, min_callers...)
variants_counts <- variants %>%
    filter(vaf_Tumor >= min_tvaf) %>%
    filter(ALT_counts_Tumor >= min_alt) %>%
    filter(cov_Tumor >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | vaf_Tumor >= 5 * vaf_Normal) %>%
    filter(ifelse(mutation_type == "SNV", n_callers >= min_callers_snv, n_callers >= min_callers_indels))

# OPTIONNAL: Add NSM info (+ wt/mut epitope sequences, immunogenicity, etc.)
if (!is.null(args$nsm_annot)) {
    # Standardize chromosome names in nsm_annot
    nsm_annot$CHROM <- standardize_chr(nsm_annot$CHROM)
    variants_counts <- nsm_annot %>% mutate(coding_consequence = "non_synonymous") %>% full_join(variants_counts, by = c("CHROM", "POS", "REF", "ALT")) %>% arrange(CHROM, POS)
    variants_counts$coding_consequence <- replace_na(variants_counts$coding_consequence, "synonymous_or_noncoding")
}

# OPTIONNAL: Use bedtoolsr to filter the mutations to keep only the ones falling in the exome bed file
if (!is.null(args$bed_exome)) {
    # Standardize chromosome names in bed_exome
    bed_exome[[1]] <- standardize_chr(bed_exome[[1]])
    
    bed_mut <- variants_counts %>%
        dplyr::select(1, 2) %>%
        mutate(end = POS) %>%
        dplyr::rename(start = POS, chrom = CHROM)
    bed_mut_exome <- bt.intersect(bed_mut, bed_exome, u = TRUE) # u=TRUE to only report one entry per mutation, if at least it is in one of the bed intervals
    bed_mut_exome <- bed_mut_exome %>% dplyr::rename(CHROM = V1, POS = V2)
    bed_mut_exome$CHROM <- standardize_chr(bed_mut_exome$CHROM)
    variants_counts <- bed_mut_exome %>%
        dplyr::select(CHROM, POS) %>%
        left_join(variants_counts, by = c("CHROM", "POS"))
}

# Get the copy number from the segment files and annotate the mutation table
# Handle case when bed_exome is not provided
if (!is.null(args$bed_exome)) {
    bed_seg_Tumor <- segs %>% dplyr::select(2:6)
    intersect_Tumor <- bt.intersect(bed_mut_exome, bed_seg_Tumor, wb = T)
    intersect_Tumor <- intersect_Tumor %>%
        dplyr::select(1, 2, 7, 8) %>%
        dplyr::rename(CHROM = V1, POS = V2, nMajor = V7, nMinor = V8)
} else {
    # Create bed_mut from variants_counts when no exome bed is provided
    bed_mut <- variants_counts %>%
        dplyr::select(1, 2) %>%
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

ccf <- ccf %>% distinct()

# Export table
write.table(ccf, file = paste0(outdir, "/", patient, "_all_mutations_CCF.tsv"), row.names = F, quote = F, sep = "\t")