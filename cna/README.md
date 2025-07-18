# CNA analysis tumor and cfDNA WES data with ASCAT 

## Download WES reference files
See https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WES

> Download the above-metioned files for both hg19 and hg38 in a single tar archive [here](https://vhio365-my.sharepoint.com/:u:/g/personal/plevy_vhio_net/ERd5t9sUYV5FnddguCqq0g8B7RRJggCKKkmqC0GbStaXeg?e=WdDOlM)

Extract `ascat_refs.tar.gz` archive into the `cna` directory.
```
tar -xzf ascat_refs.tar.gz
```

The `ascat_refs` folder should look like:
```
/ascat_refs
├── G1000_allelesAll_hg19
│   ├── G1000_alleles_hg19_chr1.txt
│   ├── G1000_alleles_hg19_chr10.txt
│   ├── G1000_alleles_hg19_chr11.txt
│   ├── G1000_alleles_hg19_chr12.txt
│   ├── G1000_alleles_hg19_chr13.txt
│   ├── G1000_alleles_hg19_chr14.txt
│   ├── G1000_alleles_hg19_chr15.txt
│   ├── G1000_alleles_hg19_chr16.txt
│   ├── G1000_alleles_hg19_chr17.txt
│   ├── G1000_alleles_hg19_chr18.txt
│   ├── G1000_alleles_hg19_chr19.txt
│   ├── G1000_alleles_hg19_chr2.txt
│   ├── G1000_alleles_hg19_chr20.txt
│   ├── G1000_alleles_hg19_chr21.txt
│   ├── G1000_alleles_hg19_chr22.txt
│   ├── G1000_alleles_hg19_chr3.txt
│   ├── G1000_alleles_hg19_chr4.txt
│   ├── G1000_alleles_hg19_chr5.txt
│   ├── G1000_alleles_hg19_chr6.txt
│   ├── G1000_alleles_hg19_chr7.txt
│   ├── G1000_alleles_hg19_chr8.txt
│   ├── G1000_alleles_hg19_chr9.txt
│   └── G1000_alleles_hg19_chrX.txt
├── G1000_allelesAll_hg38
│   ├── G1000_alleles_hg38_chr1.txt
│   ├── G1000_alleles_hg38_chr10.txt
│   ├── G1000_alleles_hg38_chr11.txt
│   ├── G1000_alleles_hg38_chr12.txt
│   ├── G1000_alleles_hg38_chr13.txt
│   ├── G1000_alleles_hg38_chr14.txt
│   ├── G1000_alleles_hg38_chr15.txt
│   ├── G1000_alleles_hg38_chr16.txt
│   ├── G1000_alleles_hg38_chr17.txt
│   ├── G1000_alleles_hg38_chr18.txt
│   ├── G1000_alleles_hg38_chr19.txt
│   ├── G1000_alleles_hg38_chr2.txt
│   ├── G1000_alleles_hg38_chr20.txt
│   ├── G1000_alleles_hg38_chr21.txt
│   ├── G1000_alleles_hg38_chr22.txt
│   ├── G1000_alleles_hg38_chr3.txt
│   ├── G1000_alleles_hg38_chr4.txt
│   ├── G1000_alleles_hg38_chr5.txt
│   ├── G1000_alleles_hg38_chr6.txt
│   ├── G1000_alleles_hg38_chr7.txt
│   ├── G1000_alleles_hg38_chr8.txt
│   ├── G1000_alleles_hg38_chr9.txt
│   └── G1000_alleles_hg38_chrX.txt
├── G1000_lociAll_hg19
│   ├── G1000_loci_hg19_chr1.txt
│   ├── G1000_loci_hg19_chr10.txt
│   ├── G1000_loci_hg19_chr11.txt
│   ├── G1000_loci_hg19_chr12.txt
│   ├── G1000_loci_hg19_chr13.txt
│   ├── G1000_loci_hg19_chr14.txt
│   ├── G1000_loci_hg19_chr15.txt
│   ├── G1000_loci_hg19_chr16.txt
│   ├── G1000_loci_hg19_chr17.txt
│   ├── G1000_loci_hg19_chr18.txt
│   ├── G1000_loci_hg19_chr19.txt
│   ├── G1000_loci_hg19_chr2.txt
│   ├── G1000_loci_hg19_chr20.txt
│   ├── G1000_loci_hg19_chr21.txt
│   ├── G1000_loci_hg19_chr22.txt
│   ├── G1000_loci_hg19_chr3.txt
│   ├── G1000_loci_hg19_chr4.txt
│   ├── G1000_loci_hg19_chr5.txt
│   ├── G1000_loci_hg19_chr6.txt
│   ├── G1000_loci_hg19_chr7.txt
│   ├── G1000_loci_hg19_chr8.txt
│   ├── G1000_loci_hg19_chr9.txt
│   └── G1000_loci_hg19_chrX.txt
├── G1000_lociAll_hg38
│   ├── G1000_loci_hg38_chr1.txt
│   ├── G1000_loci_hg38_chr10.txt
│   ├── G1000_loci_hg38_chr11.txt
│   ├── G1000_loci_hg38_chr12.txt
│   ├── G1000_loci_hg38_chr13.txt
│   ├── G1000_loci_hg38_chr14.txt
│   ├── G1000_loci_hg38_chr15.txt
│   ├── G1000_loci_hg38_chr16.txt
│   ├── G1000_loci_hg38_chr17.txt
│   ├── G1000_loci_hg38_chr18.txt
│   ├── G1000_loci_hg38_chr19.txt
│   ├── G1000_loci_hg38_chr2.txt
│   ├── G1000_loci_hg38_chr20.txt
│   ├── G1000_loci_hg38_chr21.txt
│   ├── G1000_loci_hg38_chr22.txt
│   ├── G1000_loci_hg38_chr3.txt
│   ├── G1000_loci_hg38_chr4.txt
│   ├── G1000_loci_hg38_chr5.txt
│   ├── G1000_loci_hg38_chr6.txt
│   ├── G1000_loci_hg38_chr7.txt
│   ├── G1000_loci_hg38_chr8.txt
│   ├── G1000_loci_hg38_chr9.txt
│   └── G1000_loci_hg38_chrX.txt
├── GC_G1000_hg19.txt
├── GC_G1000_hg38.txt
├── RT_G1000_hg19.txt
└── RT_G1000_hg38.txt
```

⚠️ Chromosome names in the `G1000_lociAll` files contain the `chr` prefix for `hg38` version but no prefix for `hg19` version. If your bam files were aligned to a reference containing the `chr` prefix for the `hg19` version, used this command to add the `chr` prefix to all the `G1000_lociAll_hg19` files.

```
for i in {1..22} X; do
  sed -i 's/^/chr/' G1000_loci_hg19_chr${i}.txt
done
```
or for MacOS
```
for i in {1..22} X; do
  sed -i '' 's/^/chr/' G1000_loci_hg19_chr${i}.txt
done
```

## Run ASCAT (Docker version)
📁 Run this ASCAT wrapper from within the `cna` subdirectory
```
cd cna

bash RUN_docker.sh \
	PATIENT \
	THREADS \
	SEX \
	TBAM_PATH \
	NBAM_PATH \
	OUTDIR \
	DNA \
	GENOME \
	BED_PATH \
	SKIP_NORMAL_PROCESS
```
🐳 For the **Docker** version make sure the paths (bam and bed files) are **absolute paths**. 

**Usage**: <br>
`PATIENT` : patient / sample name <br>
`THREADS` : number of CPUs to use when running ASCAT, e.g. 12 <br>
`SEX` : `XX` or `XY` (Run `bash checksex bam_file` if unknown) <br>
`TBAM_PATH` : path to tumor bam file <br>
`NBAM_PATH` : path to normal/germline bam file <br>
`OUTDIR` : desired output dir <br>
`DNA` : type of DNA source (`FrTu`, for tumor biopsy / `cfDNA`) <br>
`GENOME` : `hg19` or `hg38` <br>
`BED_PATH` : path to exome bed file <br>
`SKIP_NORMAL_PROCESS` : set to `TRUE` one if you already processed a tumor (FrTu or cfDNA) sample for this patient in a previous run, to not have to reprocess the normal bam with alleleCount, if not, set to`FALSE`. If set to `TRUE`, the `OUTDIR` has to be the same as in the previous run.

## Run ASCAT (local)
📁 Run this ASCAT wrapper from within the `cna` subdirectory

Used this version if both the ASCAT R package and alleleCount are installed locally.

```
cd cna

bash RUN.sh \
	PATIENT \
	THREADS \
	SEX \
	TBAM_PATH \
	NBAM_PATH \
	OUTDIR \
	DNA \
	GENOME \
	BED_PATH \
	SKIP_NORMAL_PROCESS
```

See **usage** above.