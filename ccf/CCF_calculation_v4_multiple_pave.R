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
  chr_vector <- as.character(chr_vector)
  chr_vector <- gsub("^chr", "", chr_vector, ignore.case = TRUE)
  return(chr_vector)
}

# Argument parser
parser <- ArgumentParser(description = "Estimate CCF for multiple Sage/Pave tumor samples (follows variant-counts --version 4)")

parser$add_argument("--patient", type = "character", required = TRUE,
                    help = "Patient ID") # patient <- "test"
parser$add_argument("--variants_counts_path", type = "character", required = TRUE,
                    help = "Path to variant read counts TSV (from variant-counts --version 4)") # variants_counts_path <- "../test/variant-counts/test_multiple_pave_variant_counts.tsv"
parser$add_argument("--config", type = "character", required = TRUE,
                    help = "Config TSV file with columns: sample_label, purity, segs_path") # config <- "../test/ccf_config.tsv"
parser$add_argument("--outdir", type = "character", required = TRUE,
                    help = "Output directory") # outdir <- "../test/ccf"

# Filters
parser$add_argument("--min_tvaf", type = "double", default = 0.03,
                    help = "Minimum tumor VAF (default: 0.03)") # min_tvaf <- 0.03
parser$add_argument("--min_alt", type = "integer", default = 4,
                    help = "Minimum ALT reads (default: 4)") # min_alt <- 4
parser$add_argument("--min_cov", type = "integer", default = 9,
                    help = "Minimum coverage (default: 9)") # min_cov <- 9

# Optional inputs
parser$add_argument("--bed_exome", type = "character", required = FALSE,
                    help = "Exome capture BED file")
parser$add_argument("--nsm_annot", type = "character", required = FALSE,
                    help = "NSM annotation table")
parser$add_argument("--neo", type = "character", required = FALSE,
                    help = "Neo output (WiGiTS)")

# Load args
args <- parser$parse_args()

patient              <- args$patient
variants_counts_path <- args$variants_counts_path
outdir               <- args$outdir
min_tvaf             <- args$min_tvaf
min_alt              <- args$min_alt
min_cov              <- args$min_cov

# Load and validate config (sample_label | purity | segs_path)
config <- fread(args$config)
required_config_cols <- c("sample_label", "purity", "segs_path")
missing_cols <- setdiff(required_config_cols, colnames(config))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in config file:", paste(missing_cols, collapse = ", ")))
}

if (!is.null(args$nsm_annot) && !is.null(args$neo)) {
  stop("Only provide one annotation file (--nsm_annot OR --neo)")
}

if (!is.null(args$bed_exome)) bed_exome <- fread(args$bed_exome)
if (!is.null(args$nsm_annot)) nsm_annot <- fread(args$nsm_annot)
if (!is.null(args$neo))       neo       <- fread(args$neo)

system(paste0("mkdir -p ", outdir))

# --------------------------------------------------------------------------- #
# Load variant counts and compute VAF / coverage from bam2R counts
# --------------------------------------------------------------------------- #
cat(sprintf("Patient: %s  |  %d samples\n", patient, nrow(config)))
variants_counts <- fread(variants_counts_path)
variants_counts$CHROM <- standardize_chr(variants_counts$CHROM)

# VAF and coverage per tumor sample
for (label in config$sample_label) {
  ref_col <- paste0("REF_counts_", label)
  alt_col <- paste0("ALT_counts_", label)
  variants_counts[[paste0("vaf_", label)]] <-
    as.numeric(variants_counts[[alt_col]]) /
    (as.numeric(variants_counts[[ref_col]]) + as.numeric(variants_counts[[alt_col]]))
  variants_counts[[paste0("cov_", label)]] <-
    as.numeric(variants_counts[[ref_col]]) + as.numeric(variants_counts[[alt_col]])
}

# VAF and coverage for Normal
variants_counts$vaf_Normal <-
  as.numeric(variants_counts$ALT_counts_Normal) /
  (as.numeric(variants_counts$REF_counts_Normal) + as.numeric(variants_counts$ALT_counts_Normal))
variants_counts$cov_Normal <-
  as.numeric(variants_counts$REF_counts_Normal) + as.numeric(variants_counts$ALT_counts_Normal)

# --------------------------------------------------------------------------- #
# Optional annotations
# --------------------------------------------------------------------------- #
if (!is.null(args$nsm_annot)) {
  nsm_annot$CHROM <- standardize_chr(nsm_annot$CHROM)
  variants_counts <- nsm_annot %>%
    mutate(coding_consequence = "non_synonymous") %>%
    full_join(variants_counts, by = c("CHROM", "POS", "REF", "ALT")) %>%
    arrange(suppressWarnings(as.numeric(CHROM)), as.numeric(POS))
  variants_counts$coding_consequence <-
    replace_na(variants_counts$coding_consequence, "synonymous_or_noncoding")
}

if (!is.null(args$neo)) {
  variants_counts <- neo %>%
    mutate(coding_consequence = "non_synonymous") %>%
    full_join(variants_counts, by = "VariantInfo") %>%
    arrange(suppressWarnings(as.numeric(CHROM)), as.numeric(POS))
  variants_counts$coding_consequence <-
    replace_na(variants_counts$coding_consequence, "synonymous_or_noncoding")
}

# --------------------------------------------------------------------------- #
# Build bed_mut for CNA intersection (optionally restricted to exome regions)
# --------------------------------------------------------------------------- #
if (!is.null(args$bed_exome)) {
  bed_exome[[1]] <- standardize_chr(bed_exome[[1]])
  bed_mut_all <- variants_counts %>%
    dplyr::select(CHROM, POS) %>%
    distinct() %>%
    mutate(end = POS) %>%
    dplyr::rename(start = POS, chrom = CHROM)
  bed_mut_exome <- bt.intersect(bed_mut_all, bed_exome, u = TRUE)
  variants_counts <- variants_counts %>%
    semi_join(bed_mut_exome, by = c("CHROM" = "V1", "POS" = "V2"))
  bed_mut <- bed_mut_exome
} else {
  bed_mut <- variants_counts %>%
    dplyr::select(CHROM, POS) %>%
    distinct() %>%
    mutate(end = POS) %>%
    dplyr::rename(start = POS, chrom = CHROM)
}

# --------------------------------------------------------------------------- #
# Per-sample: CNA intersection, purity, CCF
# --------------------------------------------------------------------------- #
for (i in seq_len(nrow(config))) {
  label      <- config$sample_label[i]
  purity_raw <- as.character(config$purity[i])
  cat(sprintf("  [%d/%d] %s  (purity = %s)\n", i, nrow(config), label, purity_raw))

  # Load segments; convert PURPLE to ASCAT format if needed
  segs_raw <- fread(config$segs_path[i])
  if ("minorAlleleCopyNumber" %in% colnames(segs_raw)) {
    if (!exists("convert_purple_to_ascat")) source("convert_purple_to_ascat.R")
    segs <- convert_purple_to_ascat(segs_raw, sample_name = paste0(patient, "_", label))
  } else {
    segs <- segs_raw
  }
  segs$chr <- standardize_chr(segs$chr)

  # Intersect with CNA segments; keep one CN value per position to avoid duplicates
  bed_seg    <- segs %>% dplyr::select(2:6)
  intersect_cn <- bt.intersect(bed_mut, bed_seg, wb = TRUE) %>%
    dplyr::select(1, 2, 7, 8) %>%
    dplyr::rename(CHROM = V1, POS = V2,
                  !!paste0("nMajor_", label) := V7,
                  !!paste0("nMinor_", label) := V8)
  intersect_cn$CHROM <- standardize_chr(intersect_cn$CHROM)
  intersect_cn <- intersect_cn %>% distinct(CHROM, POS, .keep_all = TRUE)

  variants_counts <- variants_counts %>%
    left_join(intersect_cn, by = c("CHROM", "POS"))

  # Determine purity (numeric value or estimate from VAF peak)
  if (str_to_lower(purity_raw) == "vaf") {
    vaf_col  <- paste0("vaf_", label)
    tbl_pur  <- variants_counts %>% filter(!is.na(.data[[vaf_col]]) & .data[[vaf_col]] > 0)
    peak_vaf <- density(tbl_pur[[vaf_col]])$x[which.max(density(tbl_pur[[vaf_col]])$y)]
    purity_used <- 2 * peak_vaf
    cat(sprintf("    purity estimated from VAF peak: %.3f\n", purity_used))
  } else {
    purity_used <- as.numeric(purity_raw)
  }
  variants_counts[[paste0("purity_", label)]] <- purity_used

  # CCF = (VAF / purity) * ((1 - purity) * 2 + purity * (nMajor + nMinor))
  vaf_col    <- paste0("vaf_",    label)
  nMajor_col <- paste0("nMajor_", label)
  nMinor_col <- paste0("nMinor_", label)
  variants_counts[[paste0("ccf_", label)]] <-
    (variants_counts[[vaf_col]] / purity_used) *
    ((1 - purity_used) * 2 +
       purity_used * (variants_counts[[nMajor_col]] + variants_counts[[nMinor_col]]))
}

# --------------------------------------------------------------------------- #
# Per-sample filter flags and found_in_PASS
# --------------------------------------------------------------------------- #
for (label in config$sample_label) {
  alt_col <- paste0("ALT_counts_", label)
  cov_col <- paste0("cov_",        label)
  vaf_col <- paste0("vaf_",        label)
  variants_counts[[paste0("passes_filters_", label)]] <-
    !is.na(variants_counts[[vaf_col]]) &
    variants_counts[[vaf_col]] >= min_tvaf &
    as.numeric(variants_counts[[alt_col]]) >= min_alt &
    variants_counts[[cov_col]] >= min_cov &
    variants_counts$cov_Normal >= min_cov &
    (is.na(variants_counts$vaf_Normal) | variants_counts$vaf_Normal == 0 |
       variants_counts[[vaf_col]] >= 5 * variants_counts$vaf_Normal)
}

# found_in_PASS: comma-separated sample labels passing all filters per variant
pass_matrix <- do.call(cbind, lapply(config$sample_label, function(lbl) {
  variants_counts[[paste0("passes_filters_", lbl)]]
}))
variants_counts$found_in_PASS <- apply(pass_matrix, 1, function(row) {
  paste(config$sample_label[which(row)], collapse = ",")
})

# Keep only variants passing in at least one sample
variants_counts <- variants_counts %>%
  filter(found_in_PASS != "") %>%
  distinct()

cat(sprintf(
  "Retained: %d variants  |  %d in >1 sample  |  %d in single sample\n",
  nrow(variants_counts),
  sum(str_detect(variants_counts$found_in_PASS, ",")),
  sum(!str_detect(variants_counts$found_in_PASS, ","))
))

# --------------------------------------------------------------------------- #
# Export
# --------------------------------------------------------------------------- #
out_file <- file.path(outdir, paste0(patient, "_multiple_pave_CCF.tsv"))
write.table(variants_counts, file = out_file, row.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("Done. Written to: %s\n", out_file))
