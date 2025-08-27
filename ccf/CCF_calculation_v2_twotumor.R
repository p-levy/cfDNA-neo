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
                    help = "Patient ID") #  patient <- "test"
parser$add_argument("--variants_counts_path", type = "character", required = TRUE,
                    help = "Path to variant read counts TSV") # variants_counts_path <- "../test/variant-counts/test_FrTu_cfDNA_variants_counts.tsv"
parser$add_argument("--outdir", type = "character", required = TRUE,
                    help = "Output directory") # outdir <- "../test/ccf"
parser$add_argument("--segs_Tumor_1", type = "character", required = TRUE,
                    help = "ASCAT raw segments for Tumor 1") # segs_Tumor_1 <- fread("../test/cna/test_FrTu.segments_raw.txt")
parser$add_argument("--segs_Tumor_2", type = "character", required = TRUE,
                    help = "ASCAT raw segments for Tumor 2") # segs_Tumor_2 <- fread("../test/cna/test_cfDNA.segments_raw.txt")
parser$add_argument("--purity_Tumor_1", type = "character", required = TRUE,
                    help = "Purity value for Tumor 1 (ASCAT or 'vaf')") # purity_Tumor_1 <- 0.41 # Use ASCAT purity or Write "vaf" if peak_vaf has to be used!
parser$add_argument("--purity_Tumor_2", type = "character", required = TRUE,
                    help = "Purity value for Tumor 2 (ASCAT or 'vaf')") # purity_Tumor_2 <- 0.42 # Use ASCAT purity or Write "vaf" if peak_vaf has to be used!
parser$add_argument("--sample_type_1", type = "character", required = TRUE,
                    help = "Sample type for Tumor 1 (e.g. FrTu or cfDNA, same as indicated for variant-counts step)") # sample_type_1 <- "FrTu"
parser$add_argument("--sample_type_2", type = "character", required = TRUE,
                    help = "Sample type for Tumor 1 (e.g. FrTu or cfDNA, same as indicated for variant-counts step)") # sample_type_2 <- "cfDNA"


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
segs_Tumor_2 <- fread(args$segs_Tumor_2)
segs_Tumor_1  <- fread(args$segs_Tumor_1)
purity_Tumor_1 <- args$purity_Tumor_1
purity_Tumor_2 <- args$purity_Tumor_2
sample_type_1 <- args$sample_type_1
sample_type_2 <- args$sample_type_2

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
    mutate(vaf_Tumor_1 = as.numeric(get(paste0("ALT_counts_", sample_type_1))) / (as.numeric(get(paste0("REF_counts_", sample_type_1))) + as.numeric(get(paste0("ALT_counts_", sample_type_1))))) %>%
    mutate(vaf_Tumor_2 = as.numeric(get(paste0("ALT_counts_", sample_type_2))) / (as.numeric(get(paste0("REF_counts_", sample_type_2))) + as.numeric(get(paste0("ALT_counts_", sample_type_2))))) %>%
    mutate(vaf_Normal = as.numeric(ALT_counts_Normal) / (as.numeric(REF_counts_Normal) + as.numeric(ALT_counts_Normal)))

# Add Coverage
variants <- variants %>%
    mutate(cov_Tumor_1 = as.numeric(get(paste0("REF_counts_", sample_type_1))) + as.numeric(get(paste0("ALT_counts_", sample_type_1)))) %>%
    mutate(cov_Tumor_2 = as.numeric(get(paste0("REF_counts_", sample_type_2))) + as.numeric(get(paste0("ALT_counts_", sample_type_2)))) %>%
    mutate(cov_Normal = as.numeric(REF_counts_Normal) + as.numeric(ALT_counts_Normal))

# Add n_callers
variants <- variants %>%
    rowwise() %>%
    mutate(
        n_callers.Tumor_1 = ifelse(get(paste0("set_", sample_type_1)) == "Intersection", 4, length(str_split(get(paste0("set_", sample_type_1)), "-")[[1]])),
        n_callers.Tumor_2 = ifelse(get(paste0("set_", sample_type_2)) == "Intersection", 4, length(str_split(get(paste0("set_", sample_type_2)), "-")[[1]]))
    )

# Standardize chromosome names in variants
variants$CHROM <- standardize_chr(variants$CHROM)

# Standardize chromosome names in segs (assuming chromosome is in column 2)
segs_Tumor_1$chr <- standardize_chr(segs_Tumor_1$chr)
segs_Tumor_2$chr <- standardize_chr(segs_Tumor_2$chr)

# Apply filters (on tvaf, min_cov, min_callers...)
variants_Tumor_1 <- variants %>%
    filter(vaf_Tumor_1 >= min_tvaf) %>%
    filter(get(paste0("ALT_counts_", sample_type_1)) >= min_alt) %>%
    filter(cov_Tumor_1 >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | vaf_Tumor_1 >= 5 * vaf_Normal) %>%
    filter(ifelse(mutation_type == "SNV", n_callers.Tumor_1 >= min_callers_snv, n_callers.Tumor_1 >= min_callers_indels)) %>%
    filter(vaf_Tumor_2 < min_tvaf | get(paste0("ALT_counts_", sample_type_2)) < min_alt | (ifelse(str_detect(found_in, "SNV"), is.na(n_callers.Tumor_2) | n_callers.Tumor_2 < min_callers_snv, is.na(n_callers.Tumor_2) | n_callers.Tumor_2 < min_callers_indels)) | cov_Tumor_2 < min_cov | vaf_Tumor_2 < 5 * vaf_Normal) %>%
    mutate(DNA_source = ifelse(cov_Tumor_1 >= min_cov & cov_Tumor_2 >= min_cov & cov_Normal >= min_cov, paste0(sample_type_1, "_only"), paste0(sample_type_1, "_only_min_cov_in_2")))

variants_Tumor_2 <- variants %>%
    filter(vaf_Tumor_2 >= min_tvaf) %>%
    filter(get(paste0("ALT_counts_", sample_type_2)) >= min_alt) %>%
    filter(cov_Tumor_2 >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | vaf_Tumor_2 >= 5 * vaf_Normal) %>%
    filter(ifelse(str_detect(found_in, "SNV"), n_callers.Tumor_2 >= min_callers_snv, n_callers.Tumor_2 >= min_callers_indels)) %>%
    filter(vaf_Tumor_1 < min_tvaf | get(paste0("ALT_counts_", sample_type_1)) < min_alt | (ifelse(str_detect(found_in, "SNV"), is.na(n_callers.Tumor_1) | n_callers.Tumor_1 < min_callers_snv, is.na(n_callers.Tumor_1) | n_callers.Tumor_1 < min_callers_indels)) | cov_Tumor_1 < min_cov | vaf_Tumor_1 < 5 * vaf_Normal) %>%
    mutate(DNA_source = ifelse(cov_Tumor_1 >= min_cov & cov_Tumor_2 >= min_cov & cov_Normal >= min_cov, paste0(sample_type_2, "_only"), paste0(sample_type_2, "_only_min_cov_in_2")))

variants_shared <- variants %>%
    filter(vaf_Tumor_1 >= min_tvaf & vaf_Tumor_2 >= min_tvaf) %>%
    filter(get(paste0("ALT_counts_", sample_type_1)) >= min_alt & get(paste0("ALT_counts_", sample_type_2)) >= min_alt) %>%
    filter(cov_Tumor_1 >= min_cov & cov_Tumor_2 >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | (vaf_Tumor_1 >= 5 * vaf_Normal & vaf_Tumor_2 >= 5 * vaf_Normal)) %>%
    filter(ifelse(str_detect(found_in, "SNV"), (n_callers.Tumor_1 >= min_callers_snv & n_callers.Tumor_2 >= min_callers_snv) | (CHROM == "16" & POS == 706395), n_callers.Tumor_1 >= min_callers_indels & n_callers.Tumor_2 >= min_callers_indels)) %>%
    mutate(DNA_source = "shared")

variants_counts <- rbind(variants_Tumor_1, variants_Tumor_2, variants_shared)

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

	# Get the copy number from the segment files and annotate the mutation table
	bed_seg_Tumor_1 <- segs_Tumor_1 %>% dplyr::select(2:6)
	intersect_Tumor_1 <- bt.intersect(bed_mut_exome, bed_seg_Tumor_1, wb = T)
	intersect_Tumor_1 <- intersect_Tumor_1 %>%
		dplyr::select(1, 2, 7, 8) %>%
		dplyr::rename(CHROM = V1, POS = V2, !!paste0("nMajor_", sample_type_1) := V7, !!paste0("nMinor_", sample_type_1) := V8)

	bed_seg_Tumor_2 <- segs_Tumor_2 %>% dplyr::select(2:6)
	intersect_Tumor_2 <- bt.intersect(bed_mut_exome, bed_seg_Tumor_2, wb = T)
	intersect_Tumor_2 <- intersect_Tumor_2 %>%
		dplyr::select(1, 2, 7, 8) %>%
		dplyr::rename(CHROM = V1, POS = V2, !!paste0("nMajor_", sample_type_2) := V7, !!paste0("nMinor_", sample_type_2) := V8)
} else {
	bed_mut <- variants_counts %>%
        dplyr::select(1, 2) %>%
        mutate(end = POS) %>%
        dplyr::rename(start = POS, chrom = CHROM)

	# Get the copy number from the segment files and annotate the mutation table
	bed_seg_Tumor_1 <- segs_Tumor_1 %>% dplyr::select(2:6)
	intersect_Tumor_1 <- bt.intersect(bed_mut, bed_seg_Tumor_1, wb = T)
	intersect_Tumor_1 <- intersect_Tumor_1 %>%
		dplyr::select(1, 2, 7, 8) %>%
		dplyr::rename(CHROM = V1, POS = V2, !!paste0("nMajor_", sample_type_1) := V7, !!paste0("nMinor_", sample_type_1) := V8)

	bed_seg_Tumor_2 <- segs_Tumor_2 %>% dplyr::select(2:6)
	intersect_Tumor_2 <- bt.intersect(bed_mut, bed_seg_Tumor_2, wb = T)
	intersect_Tumor_2 <- intersect_Tumor_2 %>%
		dplyr::select(1, 2, 7, 8) %>%
		dplyr::rename(CHROM = V1, POS = V2, !!paste0("nMajor_", sample_type_2) := V7, !!paste0("nMinor_", sample_type_2) := V8)
}

# Final mutation table with CN annotated
intersect_Tumor_1$CHROM <- standardize_chr(intersect_Tumor_1$CHROM)
intersect_Tumor_2$CHROM <- standardize_chr(intersect_Tumor_2$CHROM)
variants_counts <- variants_counts %>%
    left_join(intersect_Tumor_1, by = c("CHROM", "POS")) %>%
    left_join(intersect_Tumor_2, by = c("CHROM", "POS"))

# Add CCF calculation
if (str_to_lower(purity_Tumor_1) == "vaf") {
    # Calculate purity based on VAF of most abundant cluster
    table_Tumor_1 <- variants_counts %>%
    rowwise() %>%
    filter(str_detect(found_in, sample_type_1))
    # Find x value for peak density
    peak_vaf_Tumor_1 <- density(table_Tumor_1$vaf_Tumor_1)$x[which.max(density(table_Tumor_1$vaf_Tumor_1)$y)]
    peak_vaf_Tumor_1
    purity_Tumor_1 <- 2 * peak_vaf_Tumor_1
} else {
    purity_Tumor_1 <- as.numeric(purity_Tumor_1)
}

if (str_to_lower(purity_Tumor_2) == "vaf") {
    # Calculate purity based on VAF of most abundant cluster
    table_Tumor_2 <- variants_counts %>%
    rowwise() %>%
    filter(str_detect(found_in, sample_type_2))
    # Find x value for peak density
    peak_vaf_Tumor_2 <- density(table_Tumor_2$vaf_Tumor_2)$x[which.max(density(table_Tumor_2$vaf_Tumor_2)$y)]
    peak_vaf_Tumor_2
    purity_Tumor_2 <- 2 * peak_vaf_Tumor_2
} else {
    purity_Tumor_2 <- as.numeric(purity_Tumor_2)
}

# CCF = (VAF / Purity)*((purity * Tumour_TotalCN) + (Normal_TotalCN*( 1 – Purity ))))
ccf <- variants_counts %>%
    mutate(!!paste0("ccf_", sample_type_1) := (vaf_Tumor_1 / purity_Tumor_1) * ((1 - purity_Tumor_1) * 2 + purity_Tumor_1 * (get(paste0("nMajor_", sample_type_1))  + get(paste0("nMinor_", sample_type_1)))),
    !!paste0("purity_", sample_type_1) := purity_Tumor_1) %>%
    mutate(!!paste0("ccf_", sample_type_2) := (vaf_Tumor_2 / purity_Tumor_2) * ((1 - purity_Tumor_2) * 2 + purity_Tumor_2 * (get(paste0("nMajor_", sample_type_2))  + get(paste0("nMinor_", sample_type_2)))),
    !!paste0("purity_", sample_type_2) := purity_Tumor_2)

ccf <- ccf %>% distinct()

# Export table
write.table(ccf, file = paste0(outdir, "/", patient, "_", sample_type_1, "_", sample_type_2, "_all_mutations_CCF.tsv"), row.names = F, quote = F, sep = "\t")