# Variant calls to read counts
Scripts used to generate **ALT/REF read count tables** for variants called for any sample type (tumor, cfDNA and Normal DNA), using [vcfR](https://github.com/knausb/vcfR) and the `bam2R`function from the [deepSNV](https://github.com/gerstung-lab/deepSNV) R library. 

```
variant-counts --version {1|2|3} [args...]
```
See below for **version-specific usage**.

## Version 1: Merged VCF
Process a merged VCF from multiple variant callers. VCF must contain a  `set` **INFO field** (summarizing which callers called reported variants, e.g. *set=mutect-varscan-strelka* or *set=Intersect*)

### How To Run
```
variant-counts --version 1 \
	patient_name \
	vcf_path \
	sample_type \
	bam_tumor \
	bam_normal \
	min_caller_snv \
	min_caller_indel \
	outdir \
	threads
```
**Usage**: <br>
`patient_name` : patient / sample name <br>
`vcf_path` : path to vcf file (**merged from multiple variant callers**) <br>
`sample_type` : sample type <br>
`bam_tumor` : path to tumor bam file <br>
`bam_normal` : path to normal/germline bam file <br>
`min_caller_snv` : minimum number of callers required to keep a SNV <br>
`min_caller_indel` : minimum number of callers required to keep an INDEL <br>
`outdir` : output directory <br>
`threads` : number of CPUs to use when running `bam2R` for multiple variants in parallel <br>

## Version 2: Two Tumor samples
Process **two** somatic variant calls (from optionally **merged** VCF files), generated from two different tumor samples (e.g. tumor biopsy and cfDNA or two different biopsies).

### How To Run
```
variant-counts --version 2 \
	patient_name \
	vcf_path_1 \
	sample_type_1 \
	bam_tumor_1 \
	vcf_path_2 \
	sample_type_2 \
	bam_tumor_2 \
	bam_normal \
	min_caller_snv \
	min_caller_indel \
	outdir \
	threads
```
**Usage**: <br>
`patient_name` : patient / sample name <br>
`vcf_path_1` : path to first vcf file (**merged from multiple variant callers**) <br>
`sample_type_1` : sample 1 type <br>
`bam_tumor_1` : path to first tumor bam file <br>
`vcf_path_2` : path to second vcf file (**merged from multiple variant callers**) <br>
`sample_type_2` : sample 2 type <br>
`bam_tumor_2` : path to second tumor bam file <br>
`bam_normal` : path to normal/germline bam file <br>
`min_caller_snv` : minimum number of callers required to keep a SNV <br>
`min_caller_indel` : minimum number of callers required to keep an INDEL <br>
`outdir` : output directory <br>
`threads` : number of CPUs to use when running `bam2R` for multiple variants in parallel <br>

## Version 3: Single-caller VCF (Sage/Pave)

Process a somatic variant call file (VCF) obtained from the variant caller **Sage** and annotated by **Pave** (from the [nf-core/oncoanalyser](https://github.com/nf-core/oncoanalyser) pipeline). Counts are not re-computed with `bam2R` but instead extracted directly from the VCF file relevant fields (`AD` and `DP`) and variant allele frequencies from the `AF` field. 

### How To Run
```
variant-counts --version 3 \
	patient_name \
	vcf_path \
	tumor_name_vcf \
	normal_name_vcf \
	outdir \
```

**Usage**: <br>
`patient_name` : patient / sample name <br>
`vcf_path` : path to vcf file (`sage.somatic.pave.vcf.gz`) <br>
`tumor_name_vcf` : tumor sample name as in VCF file <br>
`normal_name_vcf` : normal sample name as in VCF file <br>
`outdir` : output directory <br>