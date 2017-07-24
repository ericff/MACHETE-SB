#!/bin/bash -eu

#  scratch.sh
#  
#
#  Created by Gillian Hsieh on 8/18/15.
#

SORTINGDIR=${1}
FJDir=${2}
TASK_ID=${3}

echo "SORTINGDIR: ${SORTINGDIR}"
echo "FJDir: ${FJDir}"
echo "TASK_ID: ${TASK_ID}"


StemFile=${2}/StemList.txt
echo "StemFile: ${StemFile}"
STEM=`awk 'FNR == '${TASK_ID}' {print $1}' ${StemFile}`
echo "STEM: ${STEM}"

#if [ "$(ls -A ${1}/${STEM}*)"=2 ]
#then
#    echo "sorted ${STEM}files exists, skipping alphabetizing step"
#else
for file in ${1}/*${STEM}*
do
  FILENAME=$(basename $file)
  SORTEDNAME=sorted_${FILENAME}
  if ! [[ "$FILENAME" = *sorted* ]]
  then
    echo "About to do: head -n 2 ${file} > ${SORTINGDIR}/${SORTEDNAME}"
    head -n 2 ${file} > ${SORTINGDIR}/${SORTEDNAME}
    echo "About to do: tail -n +3 ${file} | sort -k 1 >> ${SORTINGDIR}/${SORTEDNAME}"
    tail -n +3 ${file} | sort -k 1 >> ${SORTINGDIR}/${SORTEDNAME}
  fi
done;
#fi

echo "completed Sorting step for all files in ${1}"
echo "completed Sorting step for all files in ${1}" >> ${2}/MasterError.txt
