# Phylogenetic clonal reconstruction using [CONIPHER](https://github.com/McGranahanLab/CONIPHER)

## Create CONIPHER `inputTSV` from `ccf --version` step's output 

⚠️ Requires a two-tumor sample `all_mutations_CCF.tsv`file (`ccf --version` output). 

### How To Run
```
cd phylo

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
                        Sample type for Tumor (e.g. FrTu or cfDNA)
  --sample_type_2 SAMPLE_TYPE_2
                        Sample type for Tumor (e.g. FrTu or cfDNA)
  --ploidy_Tumor_1 PLOIDY_TUMOR_1
                        Ploidy for Tumor
  --ploidy_Tumor_2 PLOIDY_TUMOR_2
                        Ploidy for Tumor
```

## Run CONIPHER (🐳 Docker version)

⚠️ Dependency: requires the [`run_conipher.R`](https://raw.githubusercontent.com/McGranahanLab/CONIPHER-wrapper/refs/heads/main/src/run_conipher.R) wrapper (from the McGranahan's lab [CONIPHER-wrapper](https://github.com/McGranahanLab/CONIPHER-wrapper) repo)

```
cd phylo

./phylo-docker path/to/run_conipher.R \
	--case_id CASE_ID \
	--prefix PREFIX \
	--out_dir OUT_DIR \
	--input_tsv_loc inputTSV \
	--min_cluster_size MIN_CLUSTER_SIZE \
	--nProcs THREADS \
	--merge_clusters TRUE/FALSE

Options:
        --case_id=CHARACTER
                Tumour ID

        --prefix=CHARACTER
                Sample prefix (⚠️ has to be common to all Tumour IDs in input tsv)

        --out_dir=CHARACTER
                Working directory where output should be saved

        --input_tsv_loc=CHARACTER
                File path to input mutation table in long format

        --min_cluster_size=CHARACTER
                Minimum number of mutations in a cluster to be considered

        --nProcs=CHARACTER
                Number of cores allocated to run script in parallel

        --merge_clusters=CHARACTER
                Should similar clusters be merged if possible

        -h, --help
                Show this help message and exit

Other options (from run_conipher.R):

        --input_seg_tsv_loc=CHARACTER
                File path to input segment table used for plotting

        --subclonal_copy_correction=CHARACTER
                Should subclonal copy number correction be used

        --only_truncal_subclonal_copy_correction=CHARACTER
                Should only truncal subclonal copy number correction be used

        --pyclone_yaml_loc=CHARACTER
                Location to a template yaml file. If null package default is used

        --multiple_test_correction=CHARACTER
                Should multiple testing correction be applied for the copy number correcting mutations

        --clean_clusters=CHARACTER
                Should clusters be cleaned and merged

        --clonal_cutOff=CHARACTER
                Lower threshold CCF to be considered clonal

        --propClonal_threshold=CHARACTER
                Proportion of cluster that needs to be considered clonal to merge

        --fix_absentCCFs=CHARACTER
                Should CCF of absent mutations be set to zero

        --driver_filter=CHARACTER
                What filter to use for drivers

        --burn_in=CHARACTER
                Burn-in for DP clustering

        --seed=CHARACTER
                Seed for pyclone

        --ccf_buffer=CHARACTER
                Buffer used for CCF calculations

        --pval_cutoff=CHARACTER
                P-value for copy number calculation

        --use_boot=CHARACTER
                Should bootstrapping be used

        --correct_cpn_clusters=CHARACTER
                Should clusters driven by copy number be removed

        --adjust_noisy_clusters=CHARACTER
                Should noisy clusters be adjusted

        --adjust_noisy_clusters_prop=CHARACTER
                What is minimum proportion of mutations should be present in a region to avoid cluster adjustment

        --min_ccf=CHARACTER
                Minimum CCF to consider a mutation as present

        --multi_trees=CHARACTER
                Should alternative trees be explored

```

### Refer to CONIPHER protocol paper for details on parameters and outputs
> [**CONIPHER: a computational framework for scalable phylogenetic reconstruction with error correction**](https://doi.org/10.1038/s41596-023-00913-9) <br>
Kristiana Grigoriadis, Ariana Huebner, Abigail Bunkum, Emma Colliver, Alexander M Frankell, Mark S Hill, Kerstin Thol, Nicolai J Birkbak, Charles Swanton, Simone Zaccaria, Nicholas McGranahan <br>
***Nat Protoc** 2024 Jan;19(1):159-183. doi: 10.1038/s41596-023-00913-9*

