#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
  cat("
How to use:
arg1: patient name
arg2: path to vcf
arg3: Tumor sample name in VCF
arg4: Normal sample name in VCF
arg5: outdir
      ")
  quit()
}

# Libraries
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("tidyverse"))

# Variables
patient_name <- args[1] # test: patient_name <- "test"
vcf_path <- args[2] # test: vcf_path <- "../test/vcf/test_sage.pave.vcf.gz"
tumor_name_vcf <- args[3] # test: tumor_name_vcf <- "TEST_T"
normal_name_vcf <- args[4] # test: tumor_name_vcf <- "TEST_N"
outdir <- args[5] # test: outdir <- "../test/variant-counts"

# Process VCF
vcf <- read.vcfR(vcf_path)
vcf_tib <- vcfR2tidy(vcf, info_only = FALSE)

variants <- vcf_tib$fix %>% dplyr::select(ChromKey, CHROM, POS, REF, ALT, FILTER, IMPACT, QUAL) %>% mutate(patient = patient_name)

# Add mutation type (SNV or INDEL)
variants <- variants %>% mutate(
    mutation_type = case_when(
      nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
      TRUE ~ "INDEL"
    )
)

# Extract VAF, ALT counts and coverage (for tumor and normal separate)
vaf_table_tumor <- vcf_tib$gt %>%
    dplyr::filter(Indiv == tumor_name_vcf) %>% 
    mutate(
    vaf_Tumor = gt_AF,
    ALT_counts_Tumor = as.numeric(gsub("\\d+,", "", gt_AD)),
    cov_Tumor = gt_DP
    ) %>% 
    dplyr::select(ChromKey, POS, vaf_Tumor, ALT_counts_Tumor, cov_Tumor)

vaf_table_normal <- vcf_tib$gt %>%
    dplyr::filter(Indiv == tumor_name_vcf) %>% 
    mutate(
    vaf_Normal = gt_AF,
    ALT_counts_Normal = as.numeric(gsub("\\d+,", "", gt_AD)),
    cov_Normal = gt_DP
    ) %>% 
    dplyr::select(ChromKey, POS, vaf_Normal, ALT_counts_Normal, cov_Normal)

# Ensure vaf_table_tumor and variants have the same number of rows
    stopifnot(nrow(vaf_table_tumor) == nrow(variants))
    # Check ordering by position
    stopifnot(all(vaf_table_tumor$POS == variants$POS))
    # Drop duplicate columns before binding
    vaf_table_tumor <- vaf_table_tumor[, !(names(vaf_table_tumor) %in% names(variants))]
    # Ensure vaf_table_tumor and variants have the same number of rows
    stopifnot(nrow(vaf_table_normal) == nrow(variants))
    # Check ordering by position
    stopifnot(all(vaf_table_normal$POS == variants$POS))
    # Drop duplicate columns before binding
    vaf_table_normal <- vaf_table_normal[, !(names(vaf_table_normal) %in% names(variants))]
    # Combine column-wise (like cbind)
    variants_merge <- bind_cols(variants, vaf_table_tumor, vaf_table_normal)
	# Add VariantInfo column
	variants_merge <- variants_merge %>% mutate(
		VariantInfo = paste(CHROM, POS, REF, ALT, sep = ":")
	) %>% dplyr::select(-ChromKey)

# Ensure the structure is consistent even if no variants are found
    if (nrow(variants_merge) == 0) {
        variants_merge <- tibble(
            CHROM = character(), POS = numeric(), REF = character(), ALT = character(),
            FILTER = character(), IMPACT = character(), QUAL = numeric(), patient = character(),
            mutation_type = character(), vaf_Tumor = numeric(), ALT_counts_Tumor = numeric(), cov_Tumor = numeric(), vaf_Normal = numeric(), ALT_counts_Normal = numeric(), cov_Normal = numeric()
        )
    }

# Write the results
write.table(variants_merge, file = paste0(outdir, "/", patient_name, "_variants_counts.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
