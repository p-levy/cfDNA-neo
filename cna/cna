#!/bin/bash

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { printf "
 Usage:
 \$1=PATIENT
 \$2=THREADS
 \$3=SEX (XX or XY)
 \$4=TBAM_PATH (path to tumor bam)
 \$5=NBAM_PATH (path to normal bam)
 \$6=OUTDIR
 \$7=DNA (cfDNA or FrTu)
 \$8=GENOME (hg19 or hg38)
 \$9=BED_PATH
 \$10=SKIP_NORMAL_PROCESS (TRUE or FALSE)
 "; exit 1; }
fi

set -uex

# VARIABLES
PATIENT=${1}
THREADS=${2}
SEX=${3}
TBAM_PATH=${4}
NBAM_PATH=${5}
OUTDIR=${6}
DNA=${7}
GENOME=${8}
BED_PATH=${9}
SKIP_NORMAL_PROCESS=${10}

# ASCAT #########################################################################

mkdir -p ${OUTDIR}

Rscript ascat_${DNA}.R \
	${TBAM_PATH} \
	${NBAM_PATH} \
	${PATIENT} \
	${BED_PATH} \
	ascat_refs \
	${THREADS} \
	${SEX} \
	${GENOME} \
	${OUTDIR} \
	${SKIP_NORMAL_PROCESS}
