#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
  cat("
How to use:
arg1: patient name
arg2: path to config TSV file
arg3: path to normal/germline BAM
arg4: outdir
arg5: n threads (for parallel bam2R)

Config TSV file format (tab-separated, with header):
  vcf_path     : path to sage.somatic.pave.vcf.gz
  sample_label : short label for this sample (e.g. FrTu, cfDNA, Biopsy1)
  bam_tumor    : path to the tumor BAM file

Read counts are re-extracted from ALL BAMs (every tumor + normal) via bam2R
for every variant in the union across all VCFs — including variants not called
in a given sample, so ALT reads are reported even for uncalled sites.
FILTER status (PASS/filtered) and IMPACT annotation are taken from the VCFs.
      ")
  quit()
}

# Libraries
suppressPackageStartupMessages(library("deepSNV"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))

# Variables
patient_name <- args[1] # test: patient_name <- "test"
config_file  <- args[2] # test: config_file  <- "../test/config_samples.tsv"
bam_normal   <- args[3] # test: bam_normal   <- "../test/bam/test_Normal.bam"
outdir       <- args[4] # test: outdir       <- "../test/variant-counts"
threads      <- as.integer(args[5]) # test: threads <- 10

# Read and validate config
config <- read.table(config_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
required_cols <- c("vcf_path", "sample_label", "bam_tumor")
missing_cols  <- setdiff(required_cols, colnames(config))
if (length(missing_cols) > 0) {
  stop(paste("Missing column(s) in config file:", paste(missing_cols, collapse = ", ")))
}
if (any(duplicated(config$sample_label))) {
  stop("Duplicate sample_label values found in config file. Each sample must have a unique label.")
}

# --------------------------------------------------------------------------- #
# Helper: count reads for a given allele at a locus via bam2R
# --------------------------------------------------------------------------- #
ntCount <- function(bam, chrom, start, end, nt) {
  counts <- bam2R(bam, chrom, start, end, q = 13)
  counts[, nt] + counts[, tolower(nt)]
}

# --------------------------------------------------------------------------- #
# Function: extract variant-level info (CHROM, POS, REF, ALT, IMPACT, FILTER)
# from a Sage/Pave VCF. Does NOT extract counts — those come from bam2R.
# Returns: VariantInfo, CHROM, POS, REF, ALT, IMPACT, patient, mutation_type,
#          FILTER_{sample_label}
# --------------------------------------------------------------------------- #
extract_vcf_info <- function(vcf_path, sample_label, patient_name) {
  vcf     <- read.vcfR(vcf_path, verbose = FALSE)
  vcf_tib <- vcfR2tidy(vcf, info_only = TRUE)  # info_only=TRUE: skip GT, only need INFO + fixed fields

  fix_cols <- c("CHROM", "POS", "REF", "ALT", "FILTER")
  if ("IMPACT" %in% colnames(vcf_tib$fix)) fix_cols <- c(fix_cols, "IMPACT")

  variants <- vcf_tib$fix %>%
    dplyr::select(all_of(fix_cols)) %>%
    mutate(
      patient       = patient_name,
      mutation_type = case_when(
        nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
        nchar(REF) == nchar(ALT)           ~ "MNV",
        TRUE                               ~ "INDEL"
      ),
      VariantInfo   = paste(CHROM, POS, REF, ALT, sep = ":")
    )
  if (!"IMPACT" %in% colnames(variants)) variants$IMPACT <- NA_character_

  variants %>% dplyr::rename(!!paste0("FILTER_", sample_label) := FILTER)
}

# --------------------------------------------------------------------------- #
# Read all VCFs and build the union variant table
# --------------------------------------------------------------------------- #
cat(sprintf("Patient: %s  |  %d samples\n", patient_name, nrow(config)))

all_vcf_info <- setNames(
  lapply(seq_len(nrow(config)), function(i) {
    cat(sprintf("  [%d/%d] Reading VCF: %s\n", i, nrow(config), config$sample_label[i]))
    extract_vcf_info(config$vcf_path[i], config$sample_label[i], patient_name)
  }),
  config$sample_label
)

cat("Building variant union...\n")

# One row per unique variant; IMPACT = first non-NA across samples
base_info <- bind_rows(lapply(names(all_vcf_info), function(label) {
  all_vcf_info[[label]] %>%
    dplyr::select(VariantInfo, CHROM, POS, REF, ALT, patient, mutation_type, IMPACT)
})) %>%
  group_by(VariantInfo) %>%
  summarise(
    CHROM         = first(CHROM),
    POS           = first(POS),
    REF           = first(REF),
    ALT           = first(ALT),
    patient       = first(patient),
    mutation_type = first(mutation_type),
    IMPACT        = first(na.omit(IMPACT)),
    .groups       = "drop"
  )

# found_in: comma-separated sample labels that called each variant
found_in_df <- bind_rows(lapply(names(all_vcf_info), function(label) {
  tibble(VariantInfo = all_vcf_info[[label]]$VariantInfo, sample = label)
})) %>%
  group_by(VariantInfo) %>%
  summarise(found_in = paste(unique(sample), collapse = ","), .groups = "drop")

base_info <- base_info %>% left_join(found_in_df, by = "VariantInfo")

# FILTER_{label} per sample: PASS/filtered if called, NA if not called in that sample
for (label in names(all_vcf_info)) {
  filter_col <- paste0("FILTER_", label)
  base_info <- base_info %>%
    left_join(
      all_vcf_info[[label]] %>% dplyr::select(VariantInfo, all_of(filter_col)),
      by = "VariantInfo"
    )
}

# Remove variants with two ALT alleles at the same locus (same as v2)
base_info <- base_info %>% filter(!str_detect(ALT, "[ATCG],[ATCG]"))

cat(sprintf(
  "Union: %d variants  |  %d in >1 sample  |  %d in single sample\n",
  nrow(base_info),
  sum(str_detect(base_info$found_in, ",")),
  sum(!str_detect(base_info$found_in, ","))
))

# --------------------------------------------------------------------------- #
# Parallel bam2R: re-extract read counts from ALL BAMs for ALL union variants
# --------------------------------------------------------------------------- #
cat(sprintf("Extracting read counts via bam2R (%d variants, %d threads)...\n",
            nrow(base_info), threads))

bam_tumor_list <- config$bam_tumor
sample_labels  <- config$sample_label

cl <- parallel::makeCluster(threads)
registerDoParallel(cl)

count_rows <- foreach(
  i         = seq_len(nrow(base_info)),
  .combine  = bind_rows,
  .packages = c("deepSNV", "dplyr", "stringr"),
  .export   = c("ntCount", "bam_tumor_list", "sample_labels", "bam_normal")
) %dopar% {

  row           <- base_info[i, ]
  chrom         <- row$CHROM
  pos           <- as.numeric(row$POS)
  ref           <- row$REF
  alt           <- row$ALT
  mutation_type <- row$mutation_type

  # Returns c(REF_count, ALT_count) for a given BAM
  get_counts <- function(bam) {
    REF_c <- if (mutation_type == "SNV") {
      ntCount(bam, chrom, pos, pos, ref)
    } else {
      # INDEL or MNV: use first base of REF at the anchor position
      ntCount(bam, chrom, pos, pos, str_sub(ref, 1, 1))
    }
    ALT_c <- if (mutation_type == "SNV") {
      ntCount(bam, chrom, pos, pos, alt)
    } else if (mutation_type == "MNV") {
      # MNV (e.g. AC>GT): count first ALT base as proxy — bam2R has no multi-base columns
      ntCount(bam, chrom, pos, pos, str_sub(alt, 1, 1))
    } else if (nchar(ref) > nchar(alt)) {
      # Deletion
      ntCount(bam, chrom, pos, pos, "DEL")
    } else {
      # Insertion
      ntCount(bam, chrom, pos, pos, "INS")
    }
    c(REF_c, ALT_c)
  }

  result <- list(VariantInfo = row$VariantInfo)

  # Tumor counts per sample
  for (j in seq_along(bam_tumor_list)) {
    cts <- get_counts(bam_tumor_list[j])
    result[[paste0("REF_counts_", sample_labels[j])]] <- cts[1]
    result[[paste0("ALT_counts_", sample_labels[j])]] <- cts[2]
  }

  # Normal counts
  cts_n <- get_counts(bam_normal)
  result[["REF_counts_Normal"]] <- cts_n[1]
  result[["ALT_counts_Normal"]] <- cts_n[2]

  as.data.frame(result, stringsAsFactors = FALSE)
}

stopCluster(cl)

# --------------------------------------------------------------------------- #
# Final table and output
# --------------------------------------------------------------------------- #
variant_counts <- base_info %>% left_join(count_rows, by = "VariantInfo")

out_file <- file.path(outdir, paste0(patient_name, "_multiple_pave_variant_counts.tsv"))
write.table(variant_counts, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat(sprintf("Done. Written to: %s\n", out_file))
