#!/usr/bin/env Rscript

# Load functions
source("plots.R")

# Create test CCF compare plot
ccf_compare_plot(
	ccf_tsv = "../test/ccf/test_full_FrTu_cfDNA_all_mutations_CCF.tsv", 
	patient = "Test",
	dna_presence_criteria = "vaf", # "vaf" or "alt_reads"
	min_vaf = 0.01, # not used if dna_presence_criteria = "alt_reads"
	min_alt_reads = 1, # not used if dna_presence_criteria = "vaf"
	limit_ccf = TRUE,
	show_immunogenic = TRUE,
	label_immunogenic = TRUE
	)

# Create test ASCAT plot
pdf(file = "test_ascat_plot.pdf", width = 5, height = 2.4)
plot_ascat_allelic_segments(
	segment_file = "../test/cna/test_FrTu.segments_raw.txt",
	nmaj_color = "#7D26CD",
    nmin_color = "#00868B",
    sample_id = NULL,
    exclude_chrXY = FALSE,
    min_seg_size = 1e6,
    offset = 0.07,
    line_width = 1.5,
    cn_cap = 5)
dev.off()

# Create test CNA heatmap
set.seed(2025)
pdf("cna_heatmap.pdf", width = 8, height = 6)
cna_heatmap(
	cna_input_csv = "../test/cna/heatmap/input.csv",
	order_by = "input" # "input", "tmb" or "scna"
	)
dev.off()
