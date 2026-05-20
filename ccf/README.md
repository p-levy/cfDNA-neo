# Cancer Cell Fraction / Variant Copy Number Calculations
Scripts used to add **variant copy number** / **cancer cell fraction (CCF)** information to each variant present in the `variant_counts.tsv` files used as input. 

**Variant copy number / CCF** is estimated using the following formula:

$$
\text{CCF} = \frac{\text{VAF}}{\text{purity}} \times \Big[ (1 - \text{purity}) \times 2 + \text{purity} \times (\text{nMajorCN} + \text{nMinorCN}) \Big]
$$

`VAF`: variant allele frequency (from `variant_counts.tsv` file) <br>
`purity`: tumor purity (from ASCAT or estimated from median VAF, for samples with no detectable CNA) <br>
`nMajorCN`: major copy number (from ASCAT segments) <br>
`nMinorCN`: minor copy number (from ASCAT segments) <br>	

Optionally a variant annotation table (specified by `--nsm_annot`) can be provided, for example to add information on variant-derived **neoepitopes** or other annotations. Minimally, this annotation table **must contain** the following columns `CHROM	POS	REF	ALT` but can contain additonal columns, for example: 
```
CHROM	POS	REF	ALT	Ensembl Gene name	Gene Name	cDNA change	AA change	Wt Epitope	Mut Epitope	immunogenic
1	2426977	G	C	ENSG00000149527	PLCH2	c.G1167C	p.E389D	KLKKAASVEEGDEGQDSPGGQSRGA	KLKKAASVEEGDDGQDSPGGQSRGA	no
1	2426977	G	C	ENSG00000149527	PLCH2	c.G1800C	p.E600D	-	-	no
1	4772495	G	A	ENSG00000196581	AJAP1	c.G565A	p.E189K	ETEFIAWGPTGDEEALESNTFPGVY	ETEFIAWGPTGDKEALESNTFPGVY	no
1	89657064	C	T	ENSG00000162654	GBP4	c.G796A	p.E266K	PTNDKQYLNHMDEVPEENLERHFLM	PTNDKQYLNHMDKVPEENLERHFLM	no
1	155629882	G	T	ENSG00000163374	YY1AP1	c.C1759A	p.P587T	VNPTSFPCPLNQPLVASSVSPLIVS	VNPTSFPCPLNQTLVASSVSPLIVS	no
1	111991946	TC	T	ENSG00000116459,ENSG00000116455	ATP5F1	c.225delC	p.I75fs	NSTSRENCFYGSILLSGCRPHAIGSQ	NSTSRENCFYGSICFPAVALMR	no
1	27092960	G	A	ENSG00000117713	ARID1A	c.G1742A	p.S581N	GTMANNSAGMAASPEMMGLGDVKLT	GTMANNSAGMAANPEMMGLGDVKLT	no
1	27878527	T	G	ENSG00000126705	AHDC1	c.A100C	p.T34P	LREPKYYPGGPPTPRPLLPTRPPAS	LREPKYYPGGPPPPRPLLPTRPPAS	no
1	214170169	G	C	ENSG00000117707	PROX1	c.G291C	p.L97F	MPFPGATIISQLLKNNMNKNGGTEP	MPFPGATIISQLFKNNMNKNGGTEP	no
```
### General use:
```
./ccf --version {1|2|3} [args...]
```
>See below for **version-specific usage**.

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
  --segs SEGS           ASCAT/PURPLE raw segments
  --purity PURITY       Purity value for Tumor sample (ASCAT value or 'vaf' for samples with no CNA)
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
                        ASCAT/PURPLE raw segments for Tumor 1
  --segs_Tumor_2 SEGS_TUMOR_2
                        ASCAT/PURPLE raw segments for Tumor 2
  --purity_Tumor_1 PURITY_TUMOR_1
                        Purity value for Tumor 1 (ASCAT value or 'vaf' for samples with no CNA)
  --purity_Tumor_2 PURITY_TUMOR_2
                        Purity value for Tumor 2 (ASCAT value or 'vaf' for samples with no CNA)
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
  --segs SEGS           ASCAT/PURPLE raw segments
  --purity PURITY       Purity value for Tumor sample (ASCAT value or 'vaf' for samples with no CNA)
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