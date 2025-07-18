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

TBAM_DIR=$(dirname ${TBAM_PATH})
NBAM_DIR=$(dirname ${NBAM_PATH})
BED_DIR=$(dirname ${BED_PATH})
TBAM=$(basename ${TBAM_PATH})
NBAM=$(basename ${NBAM_PATH})
BED=$(basename ${BED_PATH})

# ASCAT #########################################################################

mkdir -p ${OUTDIR}

docker run --rm \
	-v ${TBAM_DIR}:/tbam_dir \
	-v ${NBAM_DIR}:/nbam_dir \
	-v ${BED_PATH}:/bed_dir \
	-v ${OUTDIR}:/ascat \
	-v ${PWD}:/scripts \
	-v ascat_refs:/ref \
	-w /ascat \
    plev/ascat:3.1.2 \
    Rscript /scripts/ascat_${DNA}.R \
		/tbam_dir/${TBAM} \
		/nbam_dir/${NBAM} \
		${PATIENT} \
		/bed_dir/${BED} \
		/ref \
		${THREADS} \
		${SEX} \
		${GENOME} \
		/ascat \
		${SKIP_NORMAL_PROCESS}
