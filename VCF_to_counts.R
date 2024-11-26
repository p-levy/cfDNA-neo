#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
  cat("
How to use:
arg1: path to VCF (combined_calls.vcf)
arg2: sample type
arg3: patient name
arg4: path to tumor bam
arg5: path to germline bam
arg6: outdir
arg7: n cores (for multi-thread bam2R loop)
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
vcf_path <- args[1] # test: "/Volumes/datos_lab/hla_pipeline/processed_jonatan/cfDNA_rtoledo/IMN004/IMN004/variant_calling/IMN004_cfDNA_vs_IMN004_gDNA/merged_variants/IMN004_cfDNA_vs_IMN004_gDNA_combined_calls.vcf"
sample_type <- args[2] # test: "cfDNA"
patient_name <- args[3] # test: "IMN004"
bam_tumor <- args[4] # test: "/Volumes/datos_lab/hla_pipeline/processed_jonatan/cfDNA_rtoledo/IMN004/IMN004/preprocessing/IMN004_cfDNA/recalibrated/IMN004_cfDNA.recal.bam"
bam_normal <- args[5] # test: "/Volumes/datos_lab/hla_pipeline/processed_jonatan/cfDNA_rtoledo/IMN004/IMN004/preprocessing/IMN004_gDNA/recalibrated/IMN004_gDNA.recal.bam"
outdir <- args[6]
threads <- args[7]

# Process VCF (combined calls) to extract variant infos CHROM, POS, REF, ALT, set
vfc <- read.vcfR(vcf_path)
vfc_tib <- vcfR2tidy(vfc, info_only = TRUE)

variants <- vfc_tib$fix %>% dplyr::select(CHROM, POS, REF, ALT, set) %>% mutate(sample = sample_type, patient = patient_name)
## Filter variant called by at least 2 callers for SNV and 1 caller for indels
variants <- variants %>% rowwise() %>%
  filter(ifelse(str_detect(set, "indel"), set == "Intersection" | length(str_split(set, "-")[[1]]) >= 1,
                set == "Intersection" | length(str_split(set, "-")[[1]]) >= 2))

# Add mutation type (SNV or INDEL)
variants <- variants %>% mutate(
	mutation_type = ifelse(str_detect(set, "indel"), "INDEL", "SNV")
)

ntCount <- function(bam, chrom, start, end, nt) { 
  counts <- bam2R(bam, chrom, start, end, q=13)
  pos_count <- counts[,nt]
  neg_count <- counts[,tolower(nt)]
  return(pos_count+neg_count)
}

# Remove variants with two possibilities at the same locus: too difficult to compute atm
# TO DO Fix it e.g. ifelse(str_detect(ALT, "[ATCG],[ATCG]"), paste0(ntCount(bam_FrTu, CHROM, POS, POS, gsub(",[ATCG]","",ALT)), ",", ntCount(bam_FrTu, CHROM, POS, POS, gsub("[ATCG],","",ALT)))
variants <- variants %>% filter(!str_detect(ALT, "[ATCG],[ATCG]"))

# # Extract counts for each nt at variant position from bam file and add info in table (SINGLE THREAD)
# variants_counts <- variants %>% rowwise() %>% 
#   mutate(REF_counts = ifelse(str_detect(mutation_type, "INDEL"), ntCount(bam_tumor, CHROM, POS, POS, str_sub(REF, 1, 1)), 
# ntCount(bam_tumor, CHROM, POS, POS, REF))) %>% 
#   mutate(ALT_counts = ifelse(str_detect(ALT, "[ATCG],[ATCG]"), paste0(ntCount(bam_tumor, CHROM, POS, POS, gsub(",[ATCG]","",ALT)), ",",
#                                                                ntCount(bam_tumor, CHROM, POS, POS, gsub("[ATCG],","",ALT))),
#                                    ifelse(str_detect(mutation_type, "INDEL") & nchar(REF) > nchar(ALT), ntCount(bam_tumor, CHROM, POS, POS, "DEL"),
#                                    ifelse(str_detect(mutation_type, "INDEL") & nchar(REF) < nchar(ALT), ntCount(bam_tumor, CHROM, POS, POS, "INS"),
#                                    ntCount(bam_tumor, CHROM, POS, POS, ALT))))) %>% 
#   mutate(REF_counts_Normal = ifelse(str_detect(mutation_type, "INDEL"), ntCount(bam_normal, CHROM, POS, POS, str_sub(REF, 1, 1)), 
# ntCount(bam_normal, CHROM, POS, POS, REF))) %>% 
#   mutate(ALT_counts_Normal = ifelse(str_detect(ALT, "[ATCG],[ATCG]"), paste0(ntCount(bam_normal, CHROM, POS, POS, gsub(",[ATCG]","",ALT)), ",",
#                                                                ntCount(bam_normal, CHROM, POS, POS, gsub("[ATCG],","",ALT))),
#                                    ifelse(str_detect(mutation_type, "INDEL") & nchar(REF) > nchar(ALT), ntCount(bam_normal, CHROM, POS, POS, "DEL"),
#                                    ifelse(str_detect(mutation_type, "INDEL") & nchar(REF) < nchar(ALT), ntCount(bam_normal, CHROM, POS, POS, "INS"),
#                                    ntCount(bam_normal, CHROM, POS, POS, ALT)))))

# Parallelize extraction of counts
cl <- parallel::makeCluster(as.numeric(threads)) # Create a cluster using chosen number of cores
registerDoParallel(cl)

variants_counts <- foreach(i = 1:nrow(variants), .combine = rbind, .packages = c("deepSNV", "dplyr", "stringr")) %dopar% {
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
write.table(variants_counts, file = paste0(outdir, "/", patient_name, "_variants_counts.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)