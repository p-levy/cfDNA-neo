#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
    cat("
How to use:
arg1: patient
arg2: variant read counts tsv
arg3: outdir
arg4: ASCAT raw segs cfDNA
arg5: ASCAT raw segs FrTu
arg6: BED file exome seq
arg7: Jared or Neopred pipeline Excel/TSV file
arg8: comma separated immunogenic variant positions
arg9: purity FrTu (ASCAT value or 'vaf' to use peak vaf)
arg10: purity cfDNA (ASCAT value or 'vaf' to use peak vaf)

      ")
    quit()
}


# VARIABLES
patient <- args[1] #  patient <- "TIL012"
variants_counts_path <- args[2] # variants_counts_path <- "/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/TIL012/TIL012_cfDNA_FrTu_variants_counts.tsv"
outdir <- args[3] # outdir <- "/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/TIL012/"
# Load segment copy number data (ASCAT)
segs_cfDNA <- fread(args[4]) # segs_cfDNA <- fread("/Volumes/datos_lab/ascat/TIL012/TIL012_cfDNA.segments_raw.txt")
segs_FrTu <- fread(args[5]) # segs_FrTu <- fread("/Volumes/datos_lab/ascat/TIL012/TIL012_FrTu.segments_raw.txt")
# Load exome bed file (SureSelect V6+UTR padded)
bed_exome <- fread(args[6]) # bed_exome <- fread("/Volumes/datos_lab/References/References/intervals/hg38_V6+UTR/S07604624_hs_hg38/S07604624_Padded.noheader.bed")
pipeline_excel <- args[7] # pipeline_excel <- "/Volumes/datos_lab/hla_pipeline/processed_jonatan/NEXTGEN_TIL/TIL012/hg38_v1.0.1/overlap_merge/overlap_final.txt"
immunogenic <- args[8] # comma-sep list of immunogenic positions eg: immunogenic <- "none"
purity_FrTu <- args[9] # purity_FrTu <- 0.81 # Used ASCAT purity  Write "vaf" if peak_vaf has to be used!
purity_cfDNA <- args[10] # purity_cfDNA <- 0.32 # Used ASCAT purity #2*peak_vaf_cfDNA Write "vaf" if peak_vaf has to be used!

setwd(outdir)

# Load mutation count data
variants <- fread(variants_counts_path)

# Add VAFs
variants <- variants %>%
    mutate(vaf_FrTu = as.numeric(ALT_counts_FrTu) / (as.numeric(REF_counts_FrTu) + as.numeric(ALT_counts_FrTu))) %>%
    mutate(vaf_cfDNA = as.numeric(ALT_counts_cfDNA) / (as.numeric(REF_counts_cfDNA) + as.numeric(ALT_counts_cfDNA))) %>%
    mutate(vaf_Normal = as.numeric(ALT_counts_Normal) / (as.numeric(REF_counts_Normal) + as.numeric(ALT_counts_Normal)))

# Add Coverage
variants <- variants %>%
    mutate(cov_FrTu = as.numeric(REF_counts_FrTu) + as.numeric(ALT_counts_FrTu)) %>%
    mutate(cov_cfDNA = as.numeric(REF_counts_cfDNA) + as.numeric(ALT_counts_cfDNA)) %>%
    mutate(cov_Normal = as.numeric(REF_counts_Normal) + as.numeric(ALT_counts_Normal))

# Correct n_callers column (error in vcf_to_count script!)
if (grepl("GOI|Motri", patient)) {
    variants <- variants %>%
        rowwise() %>%
        mutate(
            n_callers.FrTu = ifelse(set.FrTu == "Intersection", 4, length(str_split(set.FrTu, "-")[[1]])),
            n_callers.cfDNA = ifelse(set.cfDNA == "Intersection", 4, length(str_split(set.cfDNA, "-")[[1]]))
        )
} else {
    variants <- variants %>% dplyr::rename(n_callers.FrTu = ncallers_FrTu, n_callers.cfDNA = ncallers_cfDNA)
}

# Filter variants
# TVAF >=3%
# Variant reads =or>4 reads
# Coverage = or > 9 in all 3 tumor and normal samples
# called by 2 or more callers both SNV and INDELS
# All mutations that normal VAF = 0
# If normal has VAF, VAF >5X in tumor compared to normal

min_tvaf <- 0.03
min_alt <- 4
min_cov <- 9
min_callers_snv <- 2
min_callers_indels <- 2

variants_FrTu <- variants %>%
    filter(vaf_FrTu >= min_tvaf) %>%
    filter(ALT_counts_FrTu >= min_alt) %>%
    filter(cov_FrTu >= min_cov & cov_cfDNA >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | vaf_FrTu >= 5 * vaf_Normal) %>%
    filter(ifelse(str_detect(found_in, "SNV"), n_callers.FrTu >= min_callers_snv, n_callers.FrTu >= min_callers_indels)) %>%
    filter(vaf_cfDNA < min_tvaf | ALT_counts_cfDNA < min_alt | (ifelse(str_detect(found_in, "SNV"), is.na(n_callers.cfDNA) | n_callers.cfDNA < min_callers_snv, is.na(n_callers.cfDNA) | n_callers.cfDNA < min_callers_indels)) | cov_cfDNA < min_cov | vaf_cfDNA < 5 * vaf_Normal) %>%
    mutate(DNA_source = "FrTu_only")

variants_cfDNA <- variants %>%
    filter(vaf_cfDNA >= min_tvaf) %>%
    filter(ALT_counts_cfDNA >= min_alt) %>%
    filter(cov_FrTu >= min_cov & cov_cfDNA >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | vaf_cfDNA >= 5 * vaf_Normal) %>%
    filter(ifelse(str_detect(found_in, "SNV"), n_callers.cfDNA >= min_callers_snv, n_callers.cfDNA >= min_callers_indels)) %>%
    filter(vaf_FrTu < min_tvaf | ALT_counts_FrTu < min_alt | (ifelse(str_detect(found_in, "SNV"), is.na(n_callers.FrTu) | n_callers.FrTu < min_callers_snv, is.na(n_callers.FrTu) | n_callers.FrTu < min_callers_indels)) | cov_FrTu < min_cov | vaf_FrTu < 5 * vaf_Normal) %>%
    mutate(DNA_source = "cfDNA_only")

variants_shared <- variants %>%
    filter(vaf_FrTu >= min_tvaf & vaf_cfDNA >= min_tvaf) %>%
    filter(ALT_counts_FrTu >= min_alt & ALT_counts_cfDNA >= min_alt) %>%
    filter(cov_FrTu >= min_cov & cov_cfDNA >= min_cov & cov_Normal >= min_cov) %>%
    filter(vaf_Normal == 0 | (vaf_FrTu >= 5 * vaf_Normal & vaf_cfDNA >= 5 * vaf_Normal)) %>%
    filter(ifelse(str_detect(found_in, "SNV"), n_callers.FrTu >= min_callers_snv & n_callers.cfDNA >= min_callers_snv, n_callers.FrTu >= min_callers_indels & n_callers.cfDNA >= min_callers_indels)) %>%
    mutate(DNA_source = "shared")

variants_counts <- rbind(variants_FrTu, variants_cfDNA, variants_shared)

# Import variant analysis from Jared/Neopred and add variant annotation information columns
if (grepl("xls", pipeline_excel)) {
    variant_annot_retained <- readxl::read_excel(pipeline_excel, sheet = 1) %>%
        dplyr::select(1:42) %>%
        mutate(`Filter rationale` = NA) %>%
        mutate(filter_merge = "retained", coding_consequence = "non_synonymous")
    variant_annot_removed <- readxl::read_excel(pipeline_excel, sheet = 2) %>%
        dplyr::select(1:43) %>%
        mutate(filter_merge = "removed", coding_consequence = "non_synonymous")
    variant_annot <- rbind(variant_annot_retained, variant_annot_removed)
    variant_annot <- variant_annot %>% mutate(variant_key_short = gsub("-.*", "", `Variant key`))
    variant_annot <- separate(variant_annot, variant_key_short, c("chrom", "pos"), ":") %>%
        separate(., `Variant key`, c("chrom_pos", "nt_change"), " ") %>%
        rowwise() %>%
        mutate(variant_key_short = ifelse(str_detect(`nt_change`, ">-"), paste0(chrom, ":", as.numeric(pos) - 1),
            paste0(chrom, ":", pos)
        )) %>%
        dplyr::select(variant_key_short, everything())
} else {
    # Import variant analysis from NeoPred pipeline and add variant annotation information columns
    variant_annot_retained <- fread(pipeline_excel) %>% mutate(filter_merge = "pass", coding_consequence = "non_synonymous")
    variant_annot_removed <- fread(gsub(".txt", "_discarded.txt", pipeline_excel)) %>% mutate(filter_merge = "discarded", coding_consequence = "non_synonymous")
    variant_annot <- rbind(variant_annot_retained, variant_annot_removed)

    # Add column variant key short to variant_annot table
    variant_annot <- variant_annot %>% mutate(variant_key_short = gsub(" .*", "", `Variant_key`))
    variant_annot <- separate(variant_annot, variant_key_short, c("chrom", "pos"), ":") %>%
        separate(., `Variant_key`, c("chrom_pos", "nt_change"), " ") %>%
        rowwise() %>%
        mutate(variant_key_short = ifelse(str_detect(`nt_change`, "[ATCG]{2,}>[ATCG]"), paste0(chrom, ":", as.numeric(pos) - 1),
            paste0(chrom, ":", pos)
        )) %>%
        dplyr::select(variant_key_short, everything())
}

# Add column variant key short to variants_counts table
variants_counts <- variants_counts %>%
    rowwise() %>%
    mutate(variant_key_short = paste0(CHROM, ":", POS))

# Add annot columns to variants_counts table
variants_counts_annot <- variants_counts %>% left_join(variant_annot, by = "variant_key_short")
variants_counts_annot$coding_consequence <- replace_na(variants_counts_annot$coding_consequence, "synonymous_or_noncoding")

# Add info about immunogenic mutations
immgen_var_key_short <- trimws(strsplit(immunogenic, ";")[[1]])
variants_counts_annot <- variants_counts_annot %>% mutate(immunogenic = ifelse(variant_key_short %in% immgen_var_key_short, "yes", "no"))

# Export table with all mutations and annotations
write.table(variants_counts_annot, file = paste0(patient, "_FrTu_cfDNA_all_mutations_annot.tsv"), row.names = F, quote = F, sep = "\t")


# Use bedtoolsr to:
## 1) filter the mutations to keep only the ones falling in the exome bed file (SureSelect V6+UTR padded)
# replace chr name if needed:
if (grepl("chr", variants_counts_annot$CHROM[1])) {
    variants_counts_annot$CHROM <- gsub("chr", "", variants_counts_annot$CHROM)
	bed_exome$V1 <- gsub("chr", "", bed_exome$V1)
}
bed_mut <- variants_counts_annot %>%
    dplyr::select(1, 2) %>%
    mutate(end = POS) %>%
    dplyr::rename(start = POS, chrom = CHROM)
bed_mut_exome <- bedtoolsr::bt.intersect(bed_mut, bed_exome, u = TRUE) # u=TRUE to only report one entry per mutation, if at least it is in one of the bed intervals

## 2) get the copy number from the segment files and annotate the mutation table
bed_seg_FrTu <- segs_FrTu %>% dplyr::select(2:6)
intersect_FrTu <- bedtoolsr::bt.intersect(bed_mut_exome, bed_seg_FrTu, wb = T)
intersect_FrTu <- intersect_FrTu %>%
    dplyr::select(1, 2, 7, 8) %>%
    dplyr::rename(CHROM = V1, POS = V2, nMajor_FrTu = V7, nMinor_FrTu = V8)

bed_seg_cfDNA <- segs_cfDNA %>% dplyr::select(2:6)
intersect_cfDNA <- bedtoolsr::bt.intersect(bed_mut_exome, bed_seg_cfDNA, wb = T)
intersect_cfDNA <- intersect_cfDNA %>%
    dplyr::select(1, 2, 7, 8) %>%
    dplyr::rename(CHROM = V1, POS = V2, nMajor_cfDNA = V7, nMinor_cfDNA = V8)

# Final mutation table with CN annotated
variants_counts_cn <- variants_counts_annot %>%
    left_join(intersect_FrTu, by = c("CHROM", "POS")) %>%
    left_join(intersect_cfDNA, by = c("CHROM", "POS"))
# Filter out mutations falling in regions of homozygous loss (because suspicion of being artefacts)
# This steps also removes entries for which nMajor is NA, meaning no entry in intersect tables, so not from exome bed
variants_counts_cn <- variants_counts_cn %>% filter(nMajor_cfDNA != 0 & nMajor_FrTu != 0)

# Calculate purity based on VAF of most abundant cluster
table_FrTu <- variants_counts_cn %>%
    rowwise() %>%
    filter(str_detect(found_in, "FrTu"))
table_cfDNA <- variants_counts_cn %>%
    rowwise() %>%
    filter(str_detect(found_in, "cfDNA"))

plot_vaf_FrTu <- ggplot(table_FrTu, aes(x = vaf_FrTu)) +
    geom_density()
plot_vaf_FrTu

plot_vaf_cfDNA <- ggplot(table_cfDNA, aes(x = vaf_cfDNA)) +
    geom_density()
plot_vaf_cfDNA

# Find x value for peak density
peak_vaf_FrTu <- density(table_FrTu$vaf_FrTu)$x[which.max(density(table_FrTu$vaf_FrTu)$y)]
peak_vaf_FrTu

peak_vaf_cfDNA <- density(table_cfDNA$vaf_cfDNA)$x[which.max(density(table_cfDNA$vaf_cfDNA)$y)]
peak_vaf_cfDNA

## Plot VAF or shared mutations
table_FrTu_shared <- table_FrTu %>% filter(DNA_source == "shared")
table_cfDNA_shared <- table_cfDNA %>% filter(DNA_source == "shared")

plot_vaf_FrTu_shared <- ggplot(table_FrTu_shared, aes(x = vaf_FrTu)) +
    geom_density()
plot_vaf_FrTu_shared

plot_vaf_cfDNA_shared <- ggplot(table_cfDNA_shared, aes(x = vaf_cfDNA)) +
    geom_density()
plot_vaf_cfDNA_shared

# Find x value for peak density (of shared mutations)
peak_vaf_FrTu_shared <- density(table_FrTu_shared$vaf_FrTu)$x[which.max(density(table_FrTu_shared$vaf_FrTu)$y)]
peak_vaf_FrTu_shared

peak_vaf_cfDNA_shared <- density(table_cfDNA_shared$vaf_cfDNA)$x[which.max(density(table_cfDNA_shared$vaf_cfDNA)$y)]
peak_vaf_cfDNA_shared

# Add CCF calculation
if (str_to_lower(purity_FrTu) == "vaf") {
    purity_FrTu <- 2 * peak_vaf_FrTu
} else {
    purity_FrTu <- as.numeric(purity_FrTu)
}

if (str_to_lower(purity_cfDNA) == "vaf") {
    purity_cfDNA <- 2 * peak_vaf_cfDNA
} else {
    purity_cfDNA <- as.numeric(purity_cfDNA)
}


# CCF = (VAF / Purity)*((purity * Tumour_TotalCN) + (Normal_TotalCN*( 1 – Purity )))) (Alex Frankell formula)
ccf <- variants_counts_cn %>%
    mutate(ccf_FrTu = (vaf_FrTu / purity_FrTu) * ((1 - purity_FrTu) * 2 + purity_FrTu * (nMajor_FrTu + nMinor_FrTu))) %>%
    mutate(ccf_cfDNA = (vaf_cfDNA / purity_cfDNA) * ((purity_cfDNA * (nMajor_cfDNA + nMinor_cfDNA)) + 2 * (1 - purity_cfDNA)))

ccf <- ccf %>% distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)


# Export inputs for Phylogic NDT

# Filter coverage >= 30 in at least one sample
ccf_filt <- ccf %>%
    rowwise() %>%
    filter(ifelse(str_detect(found_in, "FrTu"), REF_counts_FrTu + ALT_counts_FrTu >= 30, REF_counts_cfDNA + ALT_counts_cfDNA >= 30))

# MAF FrTu
phylogic_input_FrTu <- ccf_filt %>%
    dplyr::select(CHROM:ALT, REF_counts_FrTu, ALT_counts_FrTu, nMajor_FrTu, nMinor_FrTu) %>%
    dplyr::rename(
        Chromosome = CHROM,
        Start_position = POS,
        Reference_Allele = REF,
        Tumor_Seq_Allele2 = ALT,
        t_ref_count = REF_counts_FrTu,
        t_alt_count = ALT_counts_FrTu
    ) %>%
    rowwise() %>%
    mutate(
        local_cn_a1 = nMinor_FrTu,
        local_cn_a2 = nMajor_FrTu
    ) %>%
    dplyr::select(-c("nMajor_FrTu", "nMinor_FrTu"))
# local_cn_a1 is nMinor and local_cn_a2 is nMajor- This is because PhylogicNDT considers that the mutation is in the a2 allele and that it cannot have 0 copy.

write.table(phylogic_input_FrTu, file = "phylogic_input_FrTu.maf", row.names = F, quote = F, sep = "\t")

# MAF cfDNA
phylogic_input_cfDNA <- ccf_filt %>%
    dplyr::select(CHROM:ALT, REF_counts_cfDNA, ALT_counts_cfDNA, nMajor_cfDNA, nMinor_cfDNA) %>%
    dplyr::rename(
        Chromosome = CHROM,
        Start_position = POS,
        Reference_Allele = REF,
        Tumor_Seq_Allele2 = ALT,
        t_ref_count = REF_counts_cfDNA,
        t_alt_count = ALT_counts_cfDNA
    ) %>%
    rowwise() %>%
    mutate(
        local_cn_a1 = nMinor_cfDNA,
        local_cn_a2 = nMajor_cfDNA
    ) %>%
    dplyr::select(-c("nMajor_cfDNA", "nMinor_cfDNA"))

# local_cn_a1 is nMajor and local_cn_a2 is nMinor, except when nMajor is 1 and nMinor is 0, then local_cn_a2 is 1 and local_cn_a1 is 0. This is because PhylogicNDT considers that the mutation is in the a2 allele and that it cannot have 0 copy.

write.table(phylogic_input_cfDNA, file = "phylogic_input_cfDNA.maf", row.names = F, quote = F, sep = "\t")

# Create MySimulation_input.sif
sif <- data.frame(
    sample_id = c("FrTu", "cfDNA"), maf_fn = c("phylogic_input_FrTu.maf", "phylogic_input_cfDNA.maf"),
    seg_fn = c("", ""), purity = c(purity_FrTu, purity_cfDNA),
    timepoint = c(1, 2)
)
write.table(sif, file = "MySimulation_input.sif", quote = F, row.names = F, sep = "\t")



# Export table
write.table(ccf, file = paste0(patient, "_FrTu_cfDNA_all_mutations_CCF.annot.tsv"), row.names = F, quote = F, sep = "\t")
