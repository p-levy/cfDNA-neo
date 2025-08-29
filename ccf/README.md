# Cancer Cell Fraction / Variant Copy Number Calculations
Scripts used to add **variant copy number** / **cancer cell fraction (ccf)** information to each variant present in the `variant_counts.tsv` files used as input. 

```
ccf --version {1|2|3} [args...]
```
See below for **version-specific usage**.

## Version 1: Merged VCF

Follows `variant-counts --version 1` run (**Single Merged VCF**)

**Reminder**: VCF must contain a  `set` **INFO field** (summarizing which callers called reported variants, e.g. *set=mutect-varscan-strelka* or *set=Intersect*)

### How To Run
```
./ccf --version 1 [-h] \
	--patient PATIENT \
	--variants_counts_path VARIANTS_COUNTS_PATH \
	--outdir OUTDIR \
	--segs SEGS \
	--purity PURITY \
	[--min_tvaf MIN_TVAF] \
	[--min_alt MIN_ALT] \
	[--min_cov MIN_COV] \
	[--min_callers_snv MIN_CALLERS_SNV] \
	[--min_callers_indels MIN_CALLERS_INDELS] \
	[--bed_exome BED_EXOME] \
	[--nsm_annot NSM_ANNOT]

options:
  -h, --help            show this help message and exit
  --patient PATIENT     Patient ID
  --variants_counts_path VARIANTS_COUNTS_PATH
                        Path to variant read counts TSV
  --outdir OUTDIR       Output directory
  --segs SEGS           ASCAT raw segments
  --purity PURITY       Purity value for Tumor sample (ASCAT or 'vaf')
  --min_tvaf MIN_TVAF   Minimum tumor VAF (default: 0.03)
  --min_alt MIN_ALT     Minimum ALT reads (default: 4)
  --min_cov MIN_COV     Minimum coverage (default: 9)
  --min_callers_snv MIN_CALLERS_SNV
                        Minimum SNV callers (default: 2)
  --min_callers_indels MIN_CALLERS_INDELS
                        Minimum INDEL callers (default: 2)
  --bed_exome BED_EXOME
                        Exome capture BED file
  --nsm_annot NSM_ANNOT
                        NSM annotation table
```

## Version 2: Two Tumor samples
Follows `variant-counts --version 2` run (**Two VCF files from two tumor samples**)

### How To Run
```
./ccf --version 2 [-h] \
	--patient PATIENT \
	--variants_counts_path VARIANTS_COUNTS_PATH \
	--outdir OUTDIR \
	--segs_Tumor_1 SEGS_TUMOR_1 \
	--segs_Tumor_2 SEGS_TUMOR_2 \
	--purity_Tumor_1 PURITY_TUMOR_1 \
	--purity_Tumor_2 PURITY_TUMOR_2 \
	--sample_type_1 SAMPLE_TYPE_1 \
	--sample_type_2 SAMPLE_TYPE_2 \
	[--min_tvaf MIN_TVAF] \
	[--min_alt MIN_ALT] \
	[--min_cov MIN_COV] \
	[--min_callers_snv MIN_CALLERS_SNV] \
	[--min_callers_indels MIN_CALLERS_INDELS] \
	[--bed_exome BED_EXOME] \
	[--nsm_annot NSM_ANNOT]

options:
  -h, --help            show this help message and exit
  --patient PATIENT     Patient ID
  --variants_counts_path VARIANTS_COUNTS_PATH
                        Path to variant read counts TSV
  --outdir OUTDIR       Output directory
  --segs_Tumor_1 SEGS_TUMOR_1
                        ASCAT raw segments for Tumor 1
  --segs_Tumor_2 SEGS_TUMOR_2
                        ASCAT raw segments for Tumor 2
  --purity_Tumor_1 PURITY_TUMOR_1
                        Purity value for Tumor 1 (ASCAT or 'vaf')
  --purity_Tumor_2 PURITY_TUMOR_2
                        Purity value for Tumor 2 (ASCAT or 'vaf')
  --sample_type_1 SAMPLE_TYPE_1
                        Sample type for Tumor 1 (e.g. FrTu or cfDNA, same as
                        indicated for variant-counts step)
  --sample_type_2 SAMPLE_TYPE_2
                        Sample type for Tumor 1 (e.g. FrTu or cfDNA, same as
                        indicated for variant-counts step)
  --min_tvaf MIN_TVAF   Minimum tumor VAF (default: 0.03)
  --min_alt MIN_ALT     Minimum ALT reads (default: 4)
  --min_cov MIN_COV     Minimum coverage (default: 9)
  --min_callers_snv MIN_CALLERS_SNV
                        Minimum SNV callers (default: 2)
  --min_callers_indels MIN_CALLERS_INDELS
                        Minimum INDEL callers (default: 2)
  --bed_exome BED_EXOME
                        Exome capture BED file
  --nsm_annot NSM_ANNOT
                        NSM annotation table
```

## Version 3: Single-caller VCF (Oncoanalyser Sage/Pave)

Follows `variant-counts --version 3` run (**Single VCF from Oncoanalyse Sage/Pave**). **Only provide Neo output OR custom NSM annotation file, not both**. 


### How To Run
```
./ccf --version 3 [-h] \
	--patient PATIENT \
	--variants_counts_path VARIANTS_COUNTS_PATH \
	--outdir OUTDIR \
	--segs SEGS \
	--purity PURITY \
	[--min_tvaf MIN_TVAF] \
	[--min_alt MIN_ALT] \
	[--min_cov MIN_COV] \
	[--bed_exome BED_EXOME] \
	[--nsm_annot NSM_ANNOT] \
	[--neo NEO]

options:
  -h, --help            show this help message and exit
  --patient PATIENT     Patient ID
  --variants_counts_path VARIANTS_COUNTS_PATH
                        Path to variant read counts TSV
  --outdir OUTDIR       Output directory
  --segs SEGS           ASCAT raw segments
  --purity PURITY       Purity value for Tumor sample (ASCAT or 'vaf')
  --min_tvaf MIN_TVAF   Minimum tumor VAF (default: 0.03)
  --min_alt MIN_ALT     Minimum ALT reads (default: 4)
  --min_cov MIN_COV     Minimum coverage (default: 9)
  --bed_exome BED_EXOME
                        Exome capture BED file
  --nsm_annot NSM_ANNOT
                        NSM annotation table
  --neo NEO
                        Neo output (WiGiTS)
```