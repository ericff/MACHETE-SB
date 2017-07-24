#!/bin/bash -eu

#  AppendNaiveRept.sh
#  
#
#  Created by Gillian Hsieh on 2/1/16.
#
### The AppendNaiveRept.sh shell calls the AppendNaiveRept.py script.  This reads in the IndelsHistogram, BadFJ and BadFJ_ver2 files, and GLM report results and outputs all the results into a single file in /FJDir/reports/AppendedReports/<STEM>_naive_report_Appended.txt
FarJuncDir=${1}
GLMReportDir=${2}
INSTALLDIR=${3}
FJGLMReportsDir=${4}
TASK_ID=${5}

if [ $# -ge 5 ]
then
OUTPUTDIR="-o ${5}"
fi

STEMFILE=${1}/StemList.txt
STEM=`awk 'FNR == '${TASK_ID}' {print $1}' ${STEMFILE}`

##ml load python/2.7.5

# Removing last argument in call of AppendNaiveRept.py so that it will just use the default directory as previously index directory was I think sent without full path -EF
# python ${INSTALLDIR}/AppendNaiveRept.py -f ${1} -g ${2} -s ${STEM} -G ${4} ${OUTPUTDIR}
python ${INSTALLDIR}/AppendNaiveRept.py -f ${1} -g ${2} -s ${STEM} -G ${4}

echo "AppendNaiveReports.sh completed for ${STEM} -- check ${1}/reports/AppendedReports/${STEM}_naive_report_Appended.txt.    This is the last step of the MACHETE."  >> ${1}/MasterError.txt
