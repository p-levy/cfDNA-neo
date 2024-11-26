#!/bin/bash

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { printf "
 Usage:
 \$1=PATIENT
 \$2=THREADS
 "; exit 1; }
fi

set -uex

# VARIABLES
PATIENT=${1}
THREADS=${2}
NEOPRED=/Volumes/datos_lab/hla_pipeline/processed_jonatan/cfDNA_rtoledo/${PATIENT}/${PATIENT}
REF=/Volumes/datos_lab/References/References
BED=S07604514_hs_hg38_Padded.bed
CFDNA_NEO=/Users/plevy/Library/CloudStorage/OneDrive-VHIO/labagros/pierre/proj/collab/andrea/cfDNA-Neo

# ASCAT #########################################################################
# docker load -i /Volumes/datos_lab/ascat/docker_ascat.tar
mkdir -p /Volumes/datos_lab/ascat/${PATIENT}
docker run --rm \
	-v ${NEOPRED}:/neopred \
	-v /Volumes/datos_lab/ascat/${PATIENT}:/ascat \
	-v ${CFDNA_NEO}:/scripts \
	-v ${REF}:/ref \
	-w /ascat \
    plevy/ascat:latest \
    Rscript /scripts/ascat_cfDNA.R \
		/neopred/preprocessing/${PATIENT}_cfDNA/recalibrated/${PATIENT}_cfDNA.recal.bam \
		/neopred/preprocessing/${PATIENT}_gDNA/recalibrated/${PATIENT}_gDNA.recal.bam \
		${PATIENT} \
		/ref/${BED} \
		/ref/ASCAT \
		${THREADS}

# VCF to COUNT ###################################################################
mkdir -p /Volumes/datos_lab/cfDNA-Neo/${PATIENT}
Rscript --vanilla $CFDNA_NEO/VCF_to_counts.R \
	${NEOPRED}/variant_calling/${PATIENT}_cfDNA_vs_${PATIENT}_gDNA/merged_variants/${PATIENT}_cfDNA_vs_${PATIENT}_gDNA_combined_calls.vcf \
	cfDNA \
	${PATIENT} \
	${NEOPRED}/preprocessing/${PATIENT}_cfDNA/recalibrated/${PATIENT}_cfDNA.recal.bam \
	${NEOPRED}/preprocessing/${PATIENT}_gDNA/recalibrated/${PATIENT}_gDNA.recal.bam \
	/Volumes/datos_lab/cfDNA-Neo/${PATIENT} \
	${THREADS}
