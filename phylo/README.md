# Phylogenetic clonal reconstruction using CONIPHER

## Create CONIPHER `inputTSV` from `ccf --version 2` step's output 

⚠️ Requires a two-tumor sample `all_mutations_CCF.tsv`file (`ccf --version 2` output). 

### How To Run
```
./phylo-input.R \
	--patient PATIENT
	--ccf CCF \
	--outdir OUTDIR \
	--sample_type_1 SAMPLE_TYPE_1 \
	--sample_type_2 SAMPLE_TYPE_2 \
	--ploidy_Tumor_1 PLOIDY_TUMOR_1 \
	--ploidy_Tumor_2 PLOIDY_TUMOR_2

options:
  -h, --help            show this help message and exit
  --patient PATIENT     Patient ID
  --ccf CCF             all_mutations_CCF.tsv
  --outdir OUTDIR       Output directory
  --sample_type_1 SAMPLE_TYPE_1
                        Sample type for Tumor 1 (e.g. FrTu or cfDNA)
  --sample_type_2 SAMPLE_TYPE_2
                        Sample type for Tumor 1 (e.g. FrTu or cfDNA)
  --ploidy_Tumor_1 PLOIDY_TUMOR_1
                        Ploidy for Tumor 1
  --ploidy_Tumor_2 PLOIDY_TUMOR_2
                        Ploidy for Tumor 2
```