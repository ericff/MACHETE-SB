#!/bin/bash -eu

#  BowtieFJAligner.batch.sh
#  Slightly changed version of BowtieAligner.batch.sh
#
#  Created by Gillian Hsieh on 11/10/15.
#

BOWTIEPARAM=${1}
BOWTIEINDEX=${2}
fastqfile=${3}
BOWTIEOUTPUTDIR=${4}
# OUTPUTFILE=${5}

#module load bowtie/2.2.4

echo "bowtie2 ${1} -x ${2} -U ${3} -S ${4}"

# Removing ${5} from call 
bowtie2 ${1} -x ${2} -U ${3} -S ${4}
