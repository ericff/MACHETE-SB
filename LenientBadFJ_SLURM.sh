#!/bin/bash -eu

#
#  
#
#  Created by Gillian Hsieh on 1/18/16.
#

FarJuncFasta=${1}
BadFJDir=${2}
SPORKDir=${3}
INSTALLDIR=${4}

#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

##ml load python/2.7.5
python ${INSTALLDIR}/SplitFastaforBadFJ.py -i ${FarJuncFasta} -l 40 -o ${BadFJDir}

# Removing former references to ${STEM}
echo "Fasta sequences split into 40bp reads - LenientBadFJ_SLURM.sh complete" >> ${SPORKDir}/MasterError.txt
echo "Unable to tell if BadFJ and BadFJ ver2 alignments complete." >> ${SPORKDir}/MasterError.txt
