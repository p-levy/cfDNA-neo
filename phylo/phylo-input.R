#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  suppressWarnings({
    library(data.table)
    library(tidyverse)
    library(argparse)
  })
})

# Argument parser
parser <- ArgumentParser(description = "Creates inputTSV for CONIPHER from all_mutations_CCF.tsv")

parser$add_argument("--patient", type = "character", required = TRUE,
                    help = "Patient ID") #  patient <- "test"
parser$add_argument("--ccf", type = "character", required = TRUE,
                    help = "all_mutations_CCF.tsv") # ccf <- fread("../test/ccf/test_FrTu_cfDNA_all_mutations_CCF.tsv")
parser$add_argument("--outdir", type = "character", required = TRUE,
                    help = "Output directory") # outdir <- "../test/ccf"
parser$add_argument("--sample_type_1", type = "character", required = TRUE,
                    help = "Sample type for Tumor 1 (e.g. FrTu or cfDNA)") # sample_type_1 <- "FrTu"
parser$add_argument("--sample_type_2", type = "character", required = TRUE,
                    help = "Sample type for Tumor 1 (e.g. FrTu or cfDNA)") # sample_type_2 <- "cfDNA"
parser$add_argument("--ploidy_Tumor_1", type = "character", required = TRUE,
                    help = "Ploidy for Tumor 1") # ploidy_Tumor_1 <- 2
parser$add_argument("--ploidy_Tumor_2", type = "character", required = TRUE,
                    help = "Ploidy for Tumor 2") # ploidy_Tumor_2 <- 2

# Load args
args <- parser$parse_args()

patient <- args$patient
ccf <- fread(args$ccf)
outdir <- args$outdir
sample_type_1 <- args$sample_type_1
sample_type_2 <- args$sample_type_2
ploidy_Tumor_1 <- args$ploidy_Tumor_1
ploidy_Tumor_2 <- args$ploidy_Tumor_2

# Create output dir if not existing already
system(paste0("mkdir -p ", outdir))

# Create TSV
muts_1 <- ccf %>% 
    mutate(
        CASE_ID = patient,
        SAMPLE = paste0(patient, "_", sample_type_1),
        ACF = get(paste0("purity_", sample_type_1)),
        PLOIDY = ploidy_Tumor_1
    ) %>% 
    dplyr::rename(
        DEPTH = cov_Tumor_1,
        CHR = CHROM,
        REF_COUNT = !!paste0("REF_counts_", sample_type_1),
        VAR_COUNT = !!paste0("ALT_counts_", sample_type_1),
        COPY_NUMBER_A = !!paste0("nMajor_", sample_type_1),
        COPY_NUMBER_B = !!paste0("nMinor_", sample_type_1)
    ) %>% 
    dplyr::select(
        CASE_ID,
        SAMPLE,
        CHR,
        POS,
        REF,
        ALT,
        REF_COUNT,
        VAR_COUNT,
        DEPTH,
        COPY_NUMBER_A,
        COPY_NUMBER_B,
        ACF,
        PLOIDY
    )

muts_2 <- ccf %>% 
    mutate(
        CASE_ID = patient,
        SAMPLE = paste0(patient, "_", sample_type_2),
        ACF = get(paste0("purity_", sample_type_2)),
        PLOIDY = ploidy_Tumor_2
    ) %>% 
    dplyr::rename(
        DEPTH = cov_Tumor_2,
        CHR = CHROM,
        REF_COUNT = !!paste0("REF_counts_", sample_type_2),
        VAR_COUNT = !!paste0("ALT_counts_", sample_type_2),
        COPY_NUMBER_A = !!paste0("nMajor_", sample_type_2),
        COPY_NUMBER_B = !!paste0("nMinor_", sample_type_2)
    ) %>% 
    dplyr::select(
        CASE_ID,
        SAMPLE,
        CHR,
        POS,
        REF,
        ALT,
        REF_COUNT,
        VAR_COUNT,
        DEPTH,
        COPY_NUMBER_A,
        COPY_NUMBER_B,
        ACF,
        PLOIDY
    )

phylo_input <- rbind(muts_1, muts_2)

write.table(phylo_input, file = paste0(outdir, "/", patient, "_", sample_type_1, "_", sample_type_2, ".conipher_input.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)