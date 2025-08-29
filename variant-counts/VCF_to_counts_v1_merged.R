#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
  cat("
How to use:
arg1: patient name
arg2: path to vcf
arg3: sample type (e.g. FrTu or cfDNA)
arg4: path to tumor bam
arg5: path to germline bam
arg6: minimum number of callers for SNVs
arg7: minimum number of callers for INDELs
arg8: outdir
arg9: n cores (for multi-thread bam2R loop)
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
patient_name <- args[1] # test: patient_name <- "TEST"
vcf_path <- args[2] # test: vcf_path <- "../test/vcf/test_FrTu.vcf.gz"
sample_type <- args[3] # test: sample_type <- "FrTu"
bam_tumor <- args[4] # test: bam_tumor <- "../test/bam/test_cfDNA.bam"
bam_normal <- args[5] # test: bam_normal <- "../test/bam/test_Normal.bam"
min_caller_snv <- args[6] # test: min_caller_snv <- 2
min_caller_indel <- args[7] # test: min_caller_indel <- 1
outdir <- args[8] # test: outdir <- "~/tmp"
threads <- args[9] # test: threads <- 10

# Process VCF (combined calls) to extract variant infos CHROM, POS, REF, ALT, set
 vcf <- read.vcfR(vcf_path)
 vcf_tib <- vcfR2tidy( vcf, info_only = TRUE)

variants <-  vcf_tib$fix %>% dplyr::select(CHROM, POS, REF, ALT, set) %>% mutate(sample = sample_type, patient = patient_name)

# Add mutation type (SNV or INDEL)
variants <- variants %>% mutate(
	mutation_type = case_when(
      nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
      TRUE ~ "INDEL"
    )
)
# Filter variants called by minimum number of callers
variants <- variants %>% rowwise() %>%
  filter(ifelse(mutation_type == "INDEL", set == "Intersection" | length(str_split(set, "-")[[1]]) >= min_caller_indel,
                set == "Intersection" | length(str_split(set, "-")[[1]]) >= min_caller_snv))


ntCount <- function(bam, chrom, start, end, nt) { 
  counts <- bam2R(bam, chrom, start, end, q=13)
  pos_count <- counts[,nt]
  neg_count <- counts[,tolower(nt)]
  return(pos_count+neg_count)
}

# Remove variants with two possibilities at the same locus: too difficult to compute atm
# TO DO Fix it e.g. ifelse(str_detect(ALT, "[ATCG],[ATCG]"), paste0(ntCount(bam_FrTu, CHROM, POS, POS, gsub(",[ATCG]","",ALT)), ",", ntCount(bam_FrTu, CHROM, POS, POS, gsub("[ATCG],","",ALT)))
variants <- variants %>% filter(!str_detect(ALT, "[ATCG],[ATCG]"))

# Extract counts for each nt at variant position from bam file and add info in table
cl <- parallel::makeCluster(as.numeric(threads)) # Create a cluster using chosen number of cores
registerDoParallel(cl)

variant_counts <- foreach(i = 1:nrow(variants), .combine = rbind, .packages = c("deepSNV", "dplyr", "stringr")) %dopar% {
  row <- variants[i, ]
  chrom <- row$CHROM
  pos <- as.numeric(row$POS)
  ref <- row$REF
  alt <- row$ALT
  mutation_type <- row$mutation_type
  
  # Tumor counts
  REF_counts_Tumor <- ifelse(str_detect(mutation_type, "INDEL"), ntCount(bam_tumor, chrom, pos, pos, str_sub(ref, 1, 1)), 
                       ntCount(bam_tumor, chrom, pos, pos, ref))
  ALT_counts_Tumor <- ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) > nchar(alt), ntCount(bam_tumor, chrom, pos, pos, "DEL"),
                       ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) < nchar(alt), ntCount(bam_tumor, chrom, pos, pos, "INS"),
                              ntCount(bam_tumor, chrom, pos, pos, alt)))
  
  # Normal counts
  REF_counts_Normal <- ifelse(str_detect(mutation_type, "INDEL"), ntCount(bam_normal, chrom, pos, pos, str_sub(ref, 1, 1)), 
                              ntCount(bam_normal, chrom, pos, pos, ref))
  ALT_counts_Normal <- ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) > nchar(alt), ntCount(bam_normal, chrom, pos, pos, "DEL"),
                               ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) < nchar(alt), ntCount(bam_normal, chrom, pos, pos, "INS"),
                                      ntCount(bam_normal, chrom, pos, pos, alt)))
  
  # Combine counts
  cbind(row, REF_counts_Tumor, ALT_counts_Tumor, REF_counts_Normal, ALT_counts_Normal)
}

stopCluster(cl) # Stop the cluster

# Write the results
write.table(variant_counts, file = paste0(outdir, "/", patient_name, "_merged_variant_counts.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)