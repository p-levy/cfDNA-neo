#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
  cat("
How to use:
arg1: path to FrTu bam
args2: path to normal bam
args3: patient name
args4: bed file
args5: path to ascat ref dir
args6: n cores to use
args7: sex ('XX' or 'XY')
args8: genome (hg19 or hg38)
args9: outdir

      ")
  quit()
}

# Libraries
suppressPackageStartupMessages(library(ASCAT))

# Variables
tumourseqfile_path <- args[1]
normalseqfile_path <- args[2]
patient <- args[3]
bedfile_path <- args[4]
ascat_ref_path <- args[5]
threads <- args[6]
sex = args[7] # 'XX' or 'XY'
genome=args[8] # hg19 or hg38
genome_number = gsub("hg", "", genome)
bedfile_path <- gsub("hg\\d{2}", genome, bedfile_path)
outdir <- args[9]
skip_normal_process=as.logical(args[9])  # converts "TRUE"/"FALSE" to logical

setwd(outdir)

# References
alleles <- paste0(ascat_ref_path, "/G1000_allelesAll_hg", genome_number, "/G1000_alleles_hg", genome_number, "_chr")
loci <- paste0(ascat_ref_path, "/G1000_lociAll_hg", genome_number, "/G1000_loci_hg", genome_number, "_chr") # Need to change chromosome name!! for i in {1..22} X; do sed -i 's/^/chr/' G1000_loci_hg", genome_number, "_chr${i}.txt; done
gc <- paste0(ascat_ref_path, "/GC_G1000_hg", genome_number, ".txt")
rt <- paste0(ascat_ref_path, "/RT_G1000_hg", genome_number, ".txt")

ascat.prepareHTS(
  tumourseqfile = tumourseqfile_path,
  normalseqfile = normalseqfile_path,
  tumourname = paste0(patient, "_FrTu"),
  normalname = paste0(patient, "_Normal"),
  allelecounter_exe = "/opt/bin/alleleCounter",
  alleles.prefix = alleles,
  loci.prefix = loci,
  gender = sex,
  genomeVersion = genome,
  nthreads = threads,
  BED_file = bedfile_path,
  chrom_names = c(1:22,'X'),
  skip_allele_counting_tumour = FALSE,
  skip_allele_counting_normal = skip_normal_process) # set to TRUE if already done

ascat.bc = ascat.loadData(Tumor_LogR_file = paste0(patient, "_FrTu_tumourLogR.txt"), Tumor_BAF_file = paste0(patient, "_FrTu_tumourBAF.txt"), Germline_LogR_file = paste0(patient, "_FrTu_normalLogR.txt"), Germline_BAF_file = paste0(patient, "_FrTu_normalBAF.txt"), gender = sex, genomeVersion = genome)
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc, replictimingfile = rt)
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')