#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

# VARIABLES
epitope_key <- fread(args[1]) # pipeline info eg: epitope_key <- fread("/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/TIL015/TIL015_epitope_variant_key.tsv")
tested_tmg <- args[2] # tested tmg epitopes tested_tmg <- "/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/TIL015/tested_epitopes.txt"
all_mutations_annot <- fread(args[3]) # all_mutations_annot <- fread("/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/TIL015/TIL015_FrTu_cfDNA_all_mutations_annot.tsv")
outdir <- args[4] # outdir <- "/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/TIL015/"
variants_counts_path <- args[5] # variants_counts_path <- "/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/TIL015/TIL015_cfDNA_FrTu_variants_counts.tsv"
patient <- args[6] # patient <- "TIL015"


tested <- readLines(tested_tmg)

# Filter epitopes + variant key table for tested epitopes
pattern <- paste0("\\b(", paste(tested, collapse = "|"), ")\\b")
epitope_key_tested <- epitope_key %>% filter(apply(., 1, function(row) any(grepl(pattern, row))))

epitope_key_tested <- epitope_key_tested %>% mutate(variant_key_short = gsub("-.*", "", `Variant key`))
# epitope_key_tested <- separate(epitope_key_tested, variant_key_short, c("chrom", "pos"), ":") %>%
#     separate(., `Variant key`, c("chrom_pos", "nt_change"), " ") %>%
#     rowwise() %>%
#     mutate(variant_key_short = ifelse(str_detect(`nt_change`, ">-"), paste0(chrom, ":", as.numeric(pos) - 1),
#         paste0(chrom, ":", pos)
#     )) %>%
#     dplyr::select(variant_key_short, everything())

# setdiff(tested, epitope_key_tested$`Mut Epitope`)

# Filter all_mutations_annot for tested variant keys
all_mutations_annot_tested <- epitope_key_tested %>%
    inner_join(all_mutations_annot, by = "variant_key_short") %>%
    distinct(variant_key_short)


# # Search tested epitope variants in vcf variants
# Load mutation count data
variants <- fread(variants_counts_path)
variants <- variants %>% rowwise() %>% mutate(variant_key_short = paste0(CHROM, ":", POS))

vcf_mut_tested <- epitope_key_tested %>%
    inner_join(variants %>% dplyr::select(variant_key_short), by = "variant_key_short") %>%
    distinct(variant_key_short)

# tested variants NOT IN unfiltered vcf variants
notfound <- setdiff(epitope_key_tested$variant_key_short, vcf_mut_tested$variant_key_short)


# Output
line1 <- paste0("Number of tested epitopes: ", length(tested))
line2 <- paste0("Number of tested variants found in variant key table: ", nrow(epitope_key_tested %>% distinct(variant_key_short)))
line3 <- paste0("Number of tested variants in filtered variants: ", nrow(all_mutations_annot_tested))
line4 <- paste0("Number of tested variants NOT IN filtered variants: ", nrow(epitope_key_tested %>% distinct(variant_key_short)) - nrow(all_mutations_annot_tested))
line5 <- paste0("Number of tested variants in unfiltered vcf variants: ", nrow(vcf_mut_tested))
line6 <- paste0("Number of tested variants NOT IN unfiltered vcf variants: ", nrow(epitope_key_tested %>% distinct(variant_key_short)) - nrow(vcf_mut_tested))
line7 <- paste0("Number of tested variants NOT IN unfiltered vcf variants: ", nrow(epitope_key_tested %>% distinct(variant_key_short)) - nrow(vcf_mut_tested))
writeLines(c(patient, line1, line2, line3, line4, line5, line6, "\n", notfound), con = paste0(outdir, "/tested_epitopes_in_filtered_and_vcf.txt"))
