#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
    cat("
How to use:
arg1: patient name
arg2: path to VCF tumor 1 
arg3: sample type tumor 1 (e.g. FrTu or cfDNA)
arg4: path to tumor bam 1
arg5: path to VCF tumor 2 
arg6: sample type tumor 2 (e.g. FrTu or cfDNA)
arg7: path to tumor bam 2
arg8: path to germline bam
arg9: minimum number of callers for SNVs
arg10: minimum number of callers for INDELs
arg11: outdir
arg12: n cores (for multi-thread bam2R loop)
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
vcf_path_1 <- args[2] # test: vcf_path_1 <- "../test/vcf/test_FrTu.vcf.gz"
sample_type_1 <- args[3] # test: sample_type_1 <-"FrTu"
bam_tumor_1 <- args[4] # test: bam_tumor_1 <- "../test/bam/test_FrTu_DNA.bam"
vcf_path_2 <- args[5] # test: vcf_path_2 <-  "../test/vcf/test_cfDNA.vcf.gz"
sample_type_2 <- args[6] # test: sample_type_2 <- "cfDNA"
bam_tumor_2 <- args[7] # test: bam_tumor_2 <- "../test/bam/test_cfDNA.bam"
bam_normal <- args[8] # test: bam_normal <- "../test/bam/test_Normal.bam"
min_caller_snv <- args[9] # test: min_caller_snv <- 2
min_caller_indel <- args[10] # test: min_caller_indel <- 1
outdir <- args[11] #  test: outdir <- "~/tmp"
threads <- args[12] # test: threads <- 10

# Process VCF 1 (combined calls) to extract variant infos CHROM, POS, REF, ALT, set
vcf_1 <- read.vcfR(vcf_path_1)
vcf_1_tib <- vcfR2tidy(vcf_1, info_only = TRUE)
variants_1 <- vcf_1_tib$fix %>%
    dplyr::select(CHROM, POS, REF, ALT, set) %>%
    rowwise() %>%
    mutate(
        sample_1 = sample_type_1, patient = patient_name,
        mutation_type = case_when(
      nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
      TRUE ~ "INDEL"
    )
    ) %>%
    filter(ifelse(mutation_type == "INDEL", set == "Intersection" | length(str_split(set, "-")[[1]]) >= min_caller_indel,
                set == "Intersection" | length(str_split(set, "-")[[1]]) >= min_caller_snv))

# Process VCF 2 (combined calls) to extract variant infos CHROM, POS, REF, ALT, set
vcf_2 <- read.vcfR(vcf_path_2)
vcf_2_tib <- vcfR2tidy(vcf_2, info_only = TRUE)
variants_2 <- vcf_2_tib$fix %>%
    dplyr::select(CHROM, POS, REF, ALT, set) %>%
    rowwise() %>%
    mutate(
        sample_2 = sample_type_2, patient = patient_name,
        mutation_type = case_when(
      nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
      TRUE ~ "INDEL"
    )
    ) %>% 
    filter(ifelse(mutation_type == "INDEL", set == "Intersection" | length(str_split(set, "-")[[1]]) >= min_caller_indel,
                set == "Intersection" | length(str_split(set, "-")[[1]]) >= min_caller_snv))


# Make joined table with all mutations (found in samples 1 and 2)
variants <- variants_1 %>% full_join(variants_2, by = c("CHROM", "POS", "REF", "ALT", "patient", "mutation_type"), suffix = c(paste0("_", sample_type_1), paste0("_", sample_type_2)))
variants <- variants %>% mutate(found_in = paste0(sample_1, ",", sample_2))
variants$found_in <- gsub(",NA$", "", variants$found_in)
variants$found_in <- gsub("^NA,", "", variants$found_in)

ntCount <- function(bam, chrom, start, end, nt) {
    counts <- bam2R(bam, chrom, start, end, q = 13)
    pos_count <- counts[, nt]
    neg_count <- counts[, tolower(nt)]
    return(pos_count + neg_count)
}

# Remove variants with two possibilities at the same locus: too difficult to compute atm
# TO DO Fix it e.g. ifelse(str_detect(ALT, "[ATCG],[ATCG]"), paste0(ntCount(bam_1, CHROM, POS, POS, gsub(",[ATCG]","",ALT)), ",", ntCount(bam_1, CHROM, POS, POS, gsub("[ATCG],","",ALT)))
variants <- variants %>% filter(!str_detect(ALT, "[ATCG],[ATCG]"))

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

    # Tumor 1 counts
    REF_counts_Tumor_1 <- ifelse(str_detect(mutation_type, "INDEL"), ntCount(bam_tumor_1, chrom, pos, pos, str_sub(ref, 1, 1)),
        ntCount(bam_tumor_1, chrom, pos, pos, ref)
    )
    ALT_counts_Tumor_1 <- ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) > nchar(alt), ntCount(bam_tumor_1, chrom, pos, pos, "DEL"),
        ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) < nchar(alt), ntCount(bam_tumor_1, chrom, pos, pos, "INS"),
            ntCount(bam_tumor_1, chrom, pos, pos, alt)
        )
    )

    # Tumor 2 counts
    REF_counts_Tumor_2 <- ifelse(str_detect(mutation_type, "INDEL"), ntCount(bam_tumor_2, chrom, pos, pos, str_sub(ref, 1, 1)),
        ntCount(bam_tumor_2, chrom, pos, pos, ref)
    )
    ALT_counts_Tumor_2 <- ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) > nchar(alt), ntCount(bam_tumor_2, chrom, pos, pos, "DEL"),
        ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) < nchar(alt), ntCount(bam_tumor_2, chrom, pos, pos, "INS"),
            ntCount(bam_tumor_2, chrom, pos, pos, alt)
        )
    )

    # Normal counts
    REF_counts_Normal <- ifelse(str_detect(mutation_type, "INDEL"), ntCount(bam_normal, chrom, pos, pos, str_sub(ref, 1, 1)),
        ntCount(bam_normal, chrom, pos, pos, ref)
    )
    ALT_counts_Normal <- ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) > nchar(alt), ntCount(bam_normal, chrom, pos, pos, "DEL"),
        ifelse(str_detect(mutation_type, "INDEL") & nchar(ref) < nchar(alt), ntCount(bam_normal, chrom, pos, pos, "INS"),
            ntCount(bam_normal, chrom, pos, pos, alt)
        )
    )

    # Combine counts
    cbind(row, REF_counts_Tumor_1, ALT_counts_Tumor_1, REF_counts_Tumor_2, ALT_counts_Tumor_2, REF_counts_Normal, ALT_counts_Normal)
}

stopCluster(cl) # Stop the cluster

#  Change column names
colnames(variants_counts) <- c(colnames(variants), paste0("REF_counts_", sample_type_1), paste0("ALT_counts_", sample_type_1), paste0("REF_counts_", sample_type_2), paste0("ALT_counts_", sample_type_2), "REF_counts_Normal", "ALT_counts_Normal")

# Write the results
write.table(variants_counts, file = paste0(outdir, "/", patient_name, "_", sample_type_1, "_", sample_type_2, "_variants_counts.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
