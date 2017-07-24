#!/bin/bash -eu
#$ -N machete_test
#$ -R y
#$ -m ae
#$ -M horence@stanford.edu
#$ -l h_vmem=10G
#$ -l h_rt=12:00:00

module load machete-nathankw/current 


CIRCPIPE_DIR=
OUTPUT_DIR=
EXONS=
REG_INDEL_INDICES=
CIRCREF=

python ${MACHETE}/run.py \
	--circpipe-dir "${CIRCPIPE_DIR}" \
	--output-dir "${OUTPUT_DIR}" \
	--hg19Exons "${EXONS}" \
	--reg-indel-indices "${REG_INDEL_INDICES}" \
	 --circref-dir "${CIRCREF}"
