import os
import shutil
import sys
import subprocess
import glob
from argparse import ArgumentParser
import pdb
import re

# NOTE THAT AS OF JULY 2016, THIS VERSION ASSUMES THAT THERE IS ONLY ON
#   PAIR OF FASTQ FILES. SMALL CHANGES COULD FIX THIS PRESUMABLY
#   BY ADDING LOOPS OVER FILES TO 
#   THE CHANGES (MARKED BELOW) MADE IN JULY 2016.

MACHETE = os.path.dirname(__file__)
# goal is to create distant paired ends using my paired end finder
# then use output to generate far junction library
# then use output to align to bowtie

# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for

def checkProcesses(popenDict):
	"""
	Function :
	Args     : popenDict - A dict. whose keys are subprocess.Popen instances. The value of a key is also a dict. that has the following
								keys: 'stdout', 'stderr', and 'cmd'.
	"""
	for popen in popenDict:
		popen.communicate() #hangs until job finishes
	for popen in popenDict:
		retcode = popen.returncode
		cmd = popenDict[popen]['cmd']
		stdout = popenDict[popen]['stdout']
		stdout.close()
		stdout_text = open(stdout.name,'r').read()
		stderr = popenDict[popen]['stderr']
		stderr.close()
		stderr_text = open(stderr.name,'r').read()
		if retcode:
			raise Exception("Command '{cmd}' failed with return code {retcode}. Log files are {stdout} and {stderr}. {stdout} contains {stdout_text}. {stderr} contains {stderr_text}.".format(cmd=cmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name,stdout_text=stdout_text,stderr_text=stderr_text))
			# raise Exception("Command '{cmd}' failed with return code {retcode}. Log files are {stdout} and {stderr}.".format(cmd=cmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

description = "Required args: --circpipe-dir, --output-dir, --hg19Exons, --reg-indel-indices, and --circref-dir."

parser = ArgumentParser(description=description)
parser.add_argument("--circpipe-dir",required=True,dest="CIRCPIPE_DIR",help="Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)")
parser.add_argument("--output-dir",required=True,dest="OUTPUT_DIR",help="Output directory for resuts. Directory path will be created recursively if it doesn't already exist. If it exists already, the directory will be deleted then created again.")
parser.add_argument("--user-bp-dist",dest="USERBPDIST",type=int,default=1000,help="Default is %(default)s.")
#parser.add_argument("REFGENOME",help="HG19 vs HG38;could upgrade to HG38.")
parser.add_argument("--num-junc-bases",dest="NUMBASESAROUNDJUNC",type=int,default=8,help="Default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70.")
parser.add_argument("--numIndels",dest="NumIndels",type=int,default=5,help="Default is %(default)s.")
parser.add_argument("--hg19Exons",required=True,dest="EXONS",help="Path to HG19Exons. Formerly called PICKLEDIR.")
parser.add_argument("--reg-indel-indices",required=True,dest="REG_INDEL_INDICES",help="Path to files with names like hg19_junctions_reg_indels_1.1.bt2l,hg19_junctions_reg_indels_2.rev.1.bt2l ... These are sometimes in a folder called IndelIndices.")
parser.add_argument("--circref-dir",required=True,dest="CIRCREF",help="Path to reference libraries output by KNIFE - directory that contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and hg19_junctions_scrambled bowtie indices.")

args = parser.parse_args()

CIRCPIPE_DIR = args.CIRCPIPE_DIR
OUTPUT_DIR = args.OUTPUT_DIR
if not os.path.isdir(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)
#else:
#	shutil.rmtree(OUTPUT_DIR)
#	os.makedirs(OUTPUT_DIR)
USERBPDIST = args.USERBPDIST
REFGENOME = "HG19" #args.REFGENOME
NUMBASESAROUNDJUNC = args.NUMBASESAROUNDJUNC
NumIndels = args.NumIndels
EXONS = args.EXONS
REG_INDEL_INDICES = args.REG_INDEL_INDICES
CIRCREF = args.CIRCREF

#end arg parsing


ORIG_DIR = os.path.join(CIRCPIPE_DIR,"orig")
UNALIGNEDDIR = os.path.join(ORIG_DIR,"unaligned")
GLM_DIR = os.path.join(CIRCPIPE_DIR,"circReads/glmReports")
DistantPEDir = os.path.join(OUTPUT_DIR,"DistantPEFiles")
if not os.path.exists(DistantPEDir):
	os.mkdir(DistantPEDir)
FASTADIR = os.path.join(OUTPUT_DIR,"fasta")
if not os.path.exists(FASTADIR):
	os.mkdir(FASTADIR)
BOWTIE_DIR = os.path.join(OUTPUT_DIR,"BowtieIndex")
if not os.path.exists(BOWTIE_DIR):
	os.mkdir(BOWTIE_DIR)
FARJUNCDIR = os.path.join(OUTPUT_DIR,"FarJunctionAlignments")
if not os.path.exists(FARJUNCDIR):
	os.mkdir(FARJUNCDIR)
SECONDFARJUNCDIR = os.path.join(OUTPUT_DIR,"FarJuncSecondary")
if not os.path.exists(SECONDFARJUNCDIR):
	os.mkdir(SECONDFARJUNCDIR)
BadFJDir = os.path.join(OUTPUT_DIR,"BadFJ")
if not os.path.exists(BadFJDir):
	os.mkdir(BadFJDir)
StemFile = os.path.join(OUTPUT_DIR,"StemList.txt")

reportsDir = os.path.join(OUTPUT_DIR,"reports")
if not os.path.exists(reportsDir):
	os.mkdir(reportsDir)
LOG_DIR = os.path.join(OUTPUT_DIR,"err_and_out")
if not os.path.exists(LOG_DIR):
	os.mkdir(LOG_DIR)

BowtieIndelsDir = os.path.join(OUTPUT_DIR,"BowtieIndels")
if not os.path.exists(BowtieIndelsDir):
	os.mkdir(BowtieIndelsDir)
FarJuncIndelsDir = os.path.join(OUTPUT_DIR,"FarJuncIndels")
if not os.path.exists(FarJuncIndelsDir):
	os.mkdir(FarJuncIndelsDir)
AlignedIndelsDir = os.path.join(SECONDFARJUNCDIR,"AlignedIndels")
if not os.path.exists(AlignedIndelsDir):
	os.mkdir(AlignedIndelsDir)
IndelsHistogramDir = os.path.join(OUTPUT_DIR,"IndelsHistogram")
if not os.path.exists(IndelsHistogramDir):
	os.mkdir(IndelsHistogramDir)
AppendedReportsDir = os.path.join(OUTPUT_DIR,"reports/AppendedReports")
if not os.path.exists(AppendedReportsDir):
	os.mkdir(AppendedReportsDir)
appendGlmDir = os.path.join(GLM_DIR,"AppendGLM")
if not os.path.isdir(appendGlmDir):
	os.mkdir(os.path.join(GLM_DIR,"AppendGLM"))
GLM_classInputDir = os.path.join(OUTPUT_DIR,"GLM_classInput")
if not os.path.exists(GLM_classInputDir):
	os.mkdir(GLM_classInputDir)

subprocess.check_call("rm -f {LOG_DIR}/*".format(LOG_DIR=LOG_DIR),shell=True)
subprocess.check_call("rm -f {LOG_DIR}/MasterError.txt".format(LOG_DIR=LOG_DIR),shell=True)

## This python script detects all the unique names for all pairs of files within a directory, eg. SRR12345, SRR123456, etc into a file called ${StemFile}
if os.path.isfile(StemFile):
	print("using existing StemList.txt")
else:
	print("generating StemList.txt from KNIFE output directory filenames")
	subprocess.check_call("python {MACHETE}/writeStemIDFiles.py -o {ORIG_DIR} -f {OUTPUT_DIR}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR),shell=True)

# counting # of times to go through the "PE matching" step - is the number of paired genome files ending in .sam /2
NUM_FILES = len(open(StemFile,"r").readlines())
print(NUM_FILES)

## if the program has been run before, there will be "sorted" reg and genome files.
## these are removed if they are present.
## All files from the original KNIFE alignments are sorted into alphabetical order because it is faster for python to identify read partners in two alphabetically sorted files than it is to find read pairs in two huge files where read IDs are scrambled.


## the shell AlphabetizeKNIFEreads.sh takes directories reg and genome, where we plan to search for mismatched paired ends, and sorts them alphabetically using the linux sort function

## sorting reg files
#j1_id
processes = {}
for index in range(1,NUM_FILES + 1):
	cmd = "{MACHETE}/AlphabetizeKNIFEreads.sh {ORIG_DIR}/reg {OUTPUT_DIR} {index}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index)
	print(cmd)
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_1sortReg.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_1sortReg.txt"),"w")
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


print("sorting reg files")


# sorting genome files
#j2_id
print("sorting genome files")
processes = {}
for index in range(1,NUM_FILES + 1):
	cmd = "{MACHETE}/AlphabetizeKNIFEreads.sh {ORIG_DIR}/genome {OUTPUT_DIR} {index}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index)
	print(cmd)
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_1sortGenome.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_1sortGenome.txt"),"w")
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
#
#
## finding mismatched paired end reads
## The shell PEfinder.sh takes the KNIFE alignment directory (OrigDir -1 ), the output directory (FJDir -2 ), the distance beyond which the user would consider alignments to be "discordant" (BP_distance -3 ), and the MACHETE installation directory (4) and calls a python script PEfinder_genomeAndReg_ENCODE.py.
## The python script identifies paired R1 and R2 from the genome and reg alignments and if the alignment location is > user defined bases apart then records them in an output file within the output directory: FarJunctionDirectory/DistantPEFiles/<STEM>_distant_pairs.txt
## If, for example, a read pair was found to be discordant, and R1= chrA:X, R2=chrB:Y, then the distant_pairs.txt file would contain the readID and chrA:M-N, chrB:P-Q where M-N is a window of 10,000 bases on each side of X and P-Q is a window of 10,000 bases on each side of Y.
## The window of 10,000 bases can be set by the user in the shell PEfinder.sh
#j3_id
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_2PEfinder.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_2PEfinder.txt"),"w")
	cmd = "{MACHETE}/PEfinder.sh {ORIG_DIR} {OUTPUT_DIR} {USERBPDIST} {MACHETE} {index}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR,USERBPDIST=USERBPDIST,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
	

print("Outputting Mismatched paired ends")
## Because there are lot of repeat locations in the _distant_pairs.txt file generated above, the DistantPE_Counter.py script is called by the shell of the same name to 1) eliminate duplicate locations.  Another early problem was that the fasta generated later was too enormous to run in a timely fashion so the distant_pairs.txt file is split by this script into 24 smaller files based on the chromosome # of the upstream partner.
## The shell DistantPE_Counter_genome_ENCODE.sh takes in the FarJunction output directory and MACHETE installation directory and outputs <FJDir>/DistantPEFiles/<STEM>/chr1,2,3,4,...,X,Y_Distant_PE_frequency.txt
#  The chrA_Distant_PE_frequency.txt files contain three columns: chrA:M-N, chrB:P-Q, and R, where R is the number of times that these two exact windows were matched together.  R could be used to cull the fasta file if it gets too large, but at this point we are still looking for junctions between exons if only one read pair aligned discordantly.
#
#j4_id
print("counting mismatch rates")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_3PEcounter.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_3PEcounter.txt"),"w")
	cmd = "{MACHETE}/DistantPE_Counter.sh {OUTPUT_DIR} {MACHETE} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

# sort mismatched PE by chromosome
## This is a simple shell script SortPairedEnds.sh to sort the chrA_Distant_PE_frequency.txt files into alphabetical order.  It takes FJDir/DistantPEFiles/chrA_Distant_PE_frequency.txt and outputs to same directory, sorted_chrA_Distant_PE_frequency.txt using the linux "sort" command.
## The reason for sorting is to increase the speed of the next step.

#j5_id
print("sorting mismatched PE files")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_4PEsort.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_4PEsort.txt"),"w")
	cmd = "{MACHETE}/SortPairedEnds.sh {OUTPUT_DIR} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

#make Junction fasta file by extracting sequence info from pickles
## Pickle files are binary files used for storage of large amounts of information in Python.  Here, GTF files have been repackaged into pickle files and store sequence, exon name, and exon location information.  Accessing pickles can be time consuming so target locations from the discordant reads have been broken down by chromosome and alphabetized.
## Loop goes through integers 1-24 where that corresponds to chromosome # (23 = X, 24 = Y). For each of those chromosomes, the shell makeJunctions.sh calls makeJunctions.py
## MakeJunctions.py takes in FJDir/DistantPEFiles/sorted__chrA_Distant_PE_frequency.txt, and outputs FJDir/fasta/<STEM>/<STEM>_chrAFarJunctions.fa.  It reads in the discordant windows, and searches the pickle file for all exons names/exon locations/ exon sequences on either the sense or antisense strands within the discordant windows.  Then it makes all pairs of possible exon junctions. All sequences are 300 base pairs long - 150 bases on each side of the breakpoint.  For exons that are fewer than 150 bases, the remainder of the 150 bases is padded with N's.  All pairs include fusions between two exons that are +/+, +/-, -/+, and -/-.  In the case that a sense and antisense exon are fused, then the sequence listed is the exact sequence that would occur if the fusion was read from the 5'->3' direction.  If the generated sequence was BLATted, the correct "strands" would appear in BLAT.  Similarly, for (-)/(-) exon pairs, if the generated sequence was BLATted, the exons would appear on the (-) strand in BLAT.
## THIS IS DIFFERENT IN KNIFE.  In KNIFE, if a (-)/(-) pair was BLATted, it would look as if it were +/+ because the KNIFE reverse complements (-) sequences.  Additionally in KNIFE, there is no way to detect inversions because -/+ and +/- fasta sequences are not generated.
print("make fusion fasta files")
for i in range(1,25):
	processes = {}
	if i == 23:
		i = "X"
	elif i == 24:
		i = "Y"
	for index in range(1,NUM_FILES + 1):
		stdout = open(os.path.join(LOG_DIR,str(i) + "_out_5makefasta.txt"),"w")
		stderr = open(os.path.join(LOG_DIR,str(i) + "_err_5makefasta.txt"),"w")
		cmd = "{MACHETE}/makeJunctions.sh {EXONS} {OUTPUT_DIR} {i} {MACHETE} {index}".format(MACHETE=MACHETE,EXONS=EXONS,OUTPUT_DIR=OUTPUT_DIR,i=i,index=index)
		print(cmd)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
	checkProcesses(processes)


##make single FJ fasta from all the fastas and then call bowtie indexer
##
## For each experiment, fasta files are generated for each chromosome separately as above.  The Bowtie2 call converts these into binary index files so the chromosome specific files must be concatenated into a single fasta file before generation of this index.
## The script linkfastafiles.sh uses linux to concatenate the <FJDir>/fasta/<STEM>/<STEM>_chr1,2,3,...,X,Y_FarJunctions.fa into a single large fasta <FJDir>/fasta/<STEM>_FarJunctions.fa.
## The second step of the linkfastafiles.sh calls Bowtie to build the Far Junctions bowtie index named <FJDir>/BowtieIndex/<STEM>/<STEM>_FJ_Index
#j6a_id this is B1
print("make FJ bowtie indices for each experiment")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_5FJIndexing.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_5FJIndexing.txt"),"w")
	print(cmd)
	cmd = "{MACHETE}/linkfastafiles.sh {OUTPUT_DIR} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,index=index)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

# Added July 5 2016
# align unaligned files to the FJ bowtie index
# This calls the shell AlignUnalignedtoFJ.  It takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.
#
#j8new_id this is B2
print("align unaligned reads to FJ index:  - check for /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam and /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq")
print("this align of unaligned reads was added in july 2016; it takes place before the original one")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_7newAlignFJ.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_7newAlignFJ.txt"),"w")
	cmd = "{MACHETE}/AlignUnalignedtoFJ.sh {OUTPUT_DIR} {ORIG_DIR} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

#
#
###make FJ naive report
## FarJuncNaiveReport.sh is a shell script that calls the python script FarJuncNaiveReport.py to generate the "Naive Reports".  Inputs include the MACHETE output directory, paths to the KNIFE alignment files, the amount a read should overlap the junction in order to be considered a "true" junctional alignment, and the MACHETE installation directory.
## see the FarJuncNaiveReport.sh for more info on FarJuncNaiveReport.py and details about how alignments are selected as "true" or "false", and how the a p value is calculated.
## The rate of true or anomaly alignments and p values are output to FJDir/reports/<STEM>_naive_report.txt.  Specific read ID's are also tracked and information on them can be found in FJDir/reports/IDs_<STEM>.txt.


#
# Added July 5 2016
#j9new_id this is C
print("make naive rpt, new as of jul 2016, earlier than previous one ")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_8newNaiveRpt.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_8NnewaiveRpt.txt"),"w")
	cmd = "{MACHETE}/FarJuncNaiveReport.sh {OUTPUT_DIR} {ORIG_DIR} {NUMBASESAROUNDJUNC} {MACHETE} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

#
#

# Added July 5 2016
# ASSUMES THERE IS ONLY ONE PAIR OF FASTQ FILES!
i=1
stdout = open(os.path.join(LOG_DIR,str(i) + "assumingonefile_out_getStem.txt"),"w")
stderr = open(os.path.join(LOG_DIR,str(i) + "assumingonefile_err_getStem.txt"),"w")
stemCmd = "awk 'FNR == '{i}' {{print $1}}' {StemFile}".format(i=i,StemFile=StemFile)
popen = subprocess.Popen(stemCmd,stdout=stdout,stderr=stderr,shell=True)
popen.communicate()	
stdout.close()
stderr.close()
retcode = popen.returncode
if retcode:
    raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

STEM_ASSUMING_ONE_FILE = open(stdout.name,'r').read().strip()

print("STEM_ASSUMING_ONE_FILE is:" + STEM_ASSUMING_ONE_FILE + "\n")

fasta_stem_dir = os.path.join(FASTADIR,STEM_ASSUMING_ONE_FILE)


## Now call parse_to_remove_FJ.py

stdout = open(os.path.join(LOG_DIR,"parse_out.txt"),"w")
stderr = open(os.path.join(LOG_DIR,"parse_err.txt"),"w")

cmd="python {MACHETE}/parse_to_remove_FJ.py --stem {STEM_ASSUMING_ONE_FILE} --outputdir {OUTPUT_DIR}".format(MACHETE=MACHETE, STEM_ASSUMING_ONE_FILE=STEM_ASSUMING_ONE_FILE, OUTPUT_DIR=OUTPUT_DIR)

print(cmd)

subprocess.call(cmd, shell=True, stdout=stdout, stderr=stderr)

stdout.close()
stderr.close()

# Added July 5 2016
##make single FJ fasta from all the fastas and then call bowtie indexer
##
## For each experiment, fasta files are generated for each chromosome separately as above.  The Bowtie2 call converts these into binary index files so the chromosome specific files must be concatenated into a single fasta file before generation of this index.
## The script linkfastafiles.sh uses linux to concatenate the <FJDir>/fasta/<STEM>/<STEM>_chr1,2,3,...,X,Y_FarJunctions.fa into a single large fasta <FJDir>/fasta/<STEM>_FarJunctions.fa.
## The second step of the linkfastafiles.sh calls Bowtie to build the Far Junctions bowtie index named <FJDir>/BowtieIndex/<STEM>/<STEM>_FJ_Index
#j6anew_id- this is B1
print("make FJ bowtie indices for each experiment")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_5secondnewFJIndexing.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_5secondnewFJIndexing.txt"),"w")
	print(cmd)
	cmd = "{MACHETE}/linkfastafiles.sh {OUTPUT_DIR} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,index=index)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

# Added July 5 2016
# align unaligned files to the FJ bowtie index
# This calls the shell AlignUnalignedtoFJ.  It takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.
#
#j8new_id this is B2
print("align unaligned reads to FJ index:  - check for /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam and /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq")
print("this align of unaligned reads was added in july 2016; it takes place before the original one")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_7secondnewAlignFJ.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_7secondnewAlignFJ.txt"),"w")
	cmd = "{MACHETE}/AlignUnalignedtoFJ.sh {OUTPUT_DIR} {ORIG_DIR} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


## If there is homology between a FarJunctions fasta sequence and the genome or transcriptome or a linear junction or circular junction, then the fusion read is less likely.  Alignments of the FarJunctions fasta sequences to the KNIFE reference indices, genome, transcriptome, linear junctions (reg), and scrambled junctions (junc) are created with two different bowtie parameters.  Bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align. These are just aligning the FJ Fasta to the bad juncs with various alignment parameters. Any junctions aligning to here will eventually be tagged as "BadFJ=1" in the final reports whereas if junctions don't align, they will receive a "BadFJ=0" in the final reports.

genomeIndex = os.path.join(CIRCREF,"hg19_genome")
transcriptomeIndex = os.path.join(CIRCREF,"hg19_transcriptome")
regIndex = os.path.join(CIRCREF,"hg19_junctions_reg")
juncIndex = os.path.join(CIRCREF,"hg19_junctions_scrambled")

# for BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:
# A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  Alignments are found in <FJDir>/BadFJ/<STEM>/<STEM>_BadFJto<ReferenceIndex>.sam.
BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 -p 4 --np 0 --rdg 50,50 --rfg 50,50"

for i in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(i) + "_out_getStem.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(i) + "_err_getStem.txt"),"w")
	stemCmd = "awk 'FNR == '{i}' {{print $1}}' {StemFile}".format(i=i,StemFile=StemFile)
	popen = subprocess.Popen(stemCmd,stdout=stdout,stderr=stderr,shell=True)
	popen.communicate()	
	stdout.close()
	stderr.close()
	retcode = popen.returncode
	if retcode:
		raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))
	STEM = open(stdout.name,'r').read().strip()
	FarJuncFasta = glob.glob(os.path.join(FASTADIR,STEM + "*FarJunctions.fa")) #should be a list of size 1
	FarJuncFasta = FarJuncFasta[0]
	BadFJStemDir =os.path.join(BadFJDir,STEM)
	if  os.path.isdir(BadFJStemDir):
		shutil.rmtree(BadFJStemDir)
	os.mkdir(BadFJStemDir)

	BadFJver2Dir=os.path.join(OUTPUT_DIR,"BadFJ_ver2",STEM)
	if os.path.isdir(BadFJver2Dir):
		shutil.rmtree(BadFJver2Dir)
	os.makedirs(BadFJver2Dir)

	r1file=os.path.join(BadFJver2Dir,STEM + "_FarJunctions_R1.fa")
	r2file=os.path.join(BadFJver2Dir,STEM + "_FarJunctions_R2.fa")
	#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM
	# is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file 
	# -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the 
	# last 40 non-N reads from the sequence.
	
	#j7_id=
	print("Identify Bad FJ's")
	processes = {}
	stdout = open(os.path.join(LOG_DIR,str(STEM) + "_out_6BadJunc.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(STEM) + "_err_6BadJunc.txt"),"w")
	cmd = "{MACHETE}/LenientBadFJ_SLURM.sh {FarJuncFasta} {BadFJver2Dir} {OUTPUT_DIR} {MACHETE}".format(MACHETE=MACHETE,FarJuncFasta=FarJuncFasta,BadFJver2Dir=BadFJver2Dir,OUTPUT_DIR=OUTPUT_DIR)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	stdout.close()
	stderr.close()
	retcode = popen.returncode
	if retcode:
		raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))
	BadFJtoGenomeFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoGenome.sam")
	if os.path.exists(BadFJtoGenomeFile):
		print("{BadFJtoGenomeFile} exists. To realign, please manually delete this file first".format(BadFJtoGenomeFile=BadFJtoGenomeFile))
	else:
		stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
		fasta = os.path.join(FASTADIR,"{STEM}_FarJunctions.fa".format(STEM=STEM))
		cmd = "{MACHETE}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {genomeIndex} {fasta} {BadFJtoGenomeFile}".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,genomeIndex=genomeIndex,fasta=fasta,BadFJtoGenomeFile=BadFJtoGenomeFile)
		print(cmd)
		print("BadFJ to genome")
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout.close()
		stderr.close()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

	BadFJtotranscriptomeFile = os.path.join(BadFJStemDir,STEM + "__BadFJtotranscriptome.sam")
	if os.path.exists(BadFJtotranscriptomeFile):
		print("{BadFJtotranscriptomeFile} exists.  To realign, please manually delete this file first".format(BadFJtotranscriptomeFile=BadFJtotranscriptomeFile))
	else:
		stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
		cmd = "{MACHETE}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {transcriptomeIndex} {fasta} {BadFJtotranscriptomeFile}".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,transcriptomeIndex=transcriptomeIndex,fasta=fasta,BadFJtotranscriptomeFile=BadFJtotranscriptomeFile)
		print(cmd)
		print("BadFJ to transcriptome")
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout.close()
		stderr.close()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

	BadFJtoRegFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoReg.sam")
	if os.path.exists(BadFJtoRegFile):
		print("{BadFJtoRegFile} exists. To realign, please manually delete this file first.".format(BadFJtoRegFile=BadFJtoRegFile))
	else:
		stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
		cmd = "{MACHETE}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {regIndex} {fasta} {BadFJtoRegFile}".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,regIndex=regIndex,fasta=fasta,BadFJtoRegFile=BadFJtoRegFile)
		print(cmd)
		print("BadFJ to reg")
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout.close()
		stderr.close()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

	
	BadFJtoJuncFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoJunc.sam")
	if os.path.exists(BadFJtoJuncFile):
		print("{BadFJtoJuncFile} exists. To realign, please manually delete this file first".format(BadFJtoJuncFile=BadFJtoJuncFile))
	else:
		stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
		print("BadFJ to junc: ")
		cmd = "{MACHETE}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {juncIndex} {fasta} {BadFJtoJuncFile}".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,juncIndex=juncIndex,fasta=fasta,BadFJtoJuncFile=BadFJtoJuncFile)
		print(cmd)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout.close()
		stderr.close()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))



	## Read gaps are disallowed in the first version of BadJuncs.  A second version of BadJuncs was created to also find genome/reg/junc/transcriptome alignments with gapped alignments.
## For BadFJ ver2 we use bowtie to align the reads1 and 2 as if they were paired end reads from the same strand.  We impose a minimum gap of 0 between the two and a maximum gap of 50,000 bases to allow up to a 50,000 base gapped alignment.
	genomeBOWTIEPARAM = "--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x {genomeIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtoGenome.sam".format(genomeIndex=genomeIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)
	transcriptomeBOWTIEPARAM = "--no-unal --no-mixed  --no-sq -p 8 -I 0 -X 50000 -f --ff -x {transcriptomeIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtotranscriptome.sam".format(transcriptomeIndex=transcriptomeIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)
	regBOWTIEPARAM = "--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x {regIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtoReg.sam".format(regIndex=regIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)
	juncBOWTIEPARAM ="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x {juncIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtoJunc.sam".format(juncIndex=juncIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)

	
	#do the bowtie alignment for each of the BadFJ Ver 2 indices above
	BadFJtoGenomeFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtoGenome.sam")
	if os.path.exists(BadFJtoGenomeFile):
		print( "{BadFJtoGenomeFile} exists. To realign, please manually delete this file first".format(BadFJtoGenomeFile=BadFJtoGenomeFile))
	else:
		stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w") 
		cmd = "{MACHETE}/BowtieAligner_BadFJv2.sh \"{genomeBOWTIEPARAM}\"".format(MACHETE=MACHETE,genomeBOWTIEPARAM=genomeBOWTIEPARAM)
		print(cmd)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes = {}
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
		checkProcesses(processes)
	
	print("BadFJ_ver2 to genome")
	BadFJtotranscriptomeFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtotranscriptome.sam")
	if os.path.exists(BadFJtotranscriptomeFile):
		print("{BadFJtotranscriptomeFile} exists. To realign, please manually delete this file first.".format(BadFJtotranscriptomeFile=BadFJtotranscriptomeFile))
	else:
		stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w")
		cmd = "{MACHETE}/BowtieAligner_BadFJv2.sh \"{transcriptomeBOWTIEPARAM}\"".format(MACHETE=MACHETE,transcriptomeBOWTIEPARAM=transcriptomeBOWTIEPARAM)
		print(cmd)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes = {}
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
		checkProcesses(processes)
	print("BadFJ_ver2 to transcriptome")

	BadFJtoRegFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtoReg.sam")
	if os.path.exists(BadFJtoRegFile):
		print("{BadFJtoRegFile} exists. To realign, please manually delete this file first.".format(BadFJtoRegFile=BadFJtoRegFile))
	else:
		stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w")
		cmd = "{MACHETE}/BowtieAligner_BadFJv2.sh \"{regBOWTIEPARAM}\"".format(MACHETE=MACHETE,regBOWTIEPARAM=regBOWTIEPARAM)
		print(cmd)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes = {}
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
		checkProcesses(processes)
	print("BadFJ_ver2 to reg")

	BadFJtoJuncFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtoJunc.sam")	
	if os.path.exists(BadFJtoJuncFile):
		print("{BadFJtoJuncFile} exists. To realign, please manually delete this file first.".format(BadFJtoJuncFile=BadFJtoJuncFile))
	else:
		stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
		stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w")
		cmd = "{MACHETE}/BowtieAligner_BadFJv2.sh \"{juncBOWTIEPARAM}\"".format(MACHETE=MACHETE,juncBOWTIEPARAM=juncBOWTIEPARAM)
		print(cmd)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes = {}
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
		checkProcesses(processes)
	print("BadFJ_ver2 to junc")



# align unaligned files to the FJ bowtie index
# This calls the shell AlignUnalignedtoFJ.  It takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.
#
#j8_id
print("align unaligned reads to FJ index:  - check for /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam and /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_7AlignFJ.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_7AlignFJ.txt"),"w")
	cmd = "{MACHETE}/AlignUnalignedtoFJ.sh {OUTPUT_DIR} {ORIG_DIR} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


#
#
###make FJ naive report
## FarJuncNaiveReport.sh is a shell script that calls the python script FarJuncNaiveReport.py to generate the "Naive Reports".  Inputs include the MACHETE output directory, paths to the KNIFE alignment files, the amount a read should overlap the junction in order to be considered a "true" junctional alignment, and the MACHETE installation directory.
## see the FarJuncNaiveReport.sh for more info on FarJuncNaiveReport.py and details about how alignments are selected as "true" or "false", and how the a p value is calculated.
## The rate of true or anomaly alignments and p values are output to FJDir/reports/<STEM>_naive_report.txt.  Specific read ID's are also tracked and information on them can be found in FJDir/reports/IDs_<STEM>.txt.


#
#j9_id
print("make naive rpt")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_8NaiveRpt.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_8NaiveRpt.txt"),"w")
	cmd = "{MACHETE}/FarJuncNaiveReport.sh {OUTPUT_DIR} {ORIG_DIR} {NUMBASESAROUNDJUNC} {MACHETE} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

##
####### GLM #######
##Make Class input files for GLM
##
###
#j15a_id
print("make FJ class input files")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_15FJforGLM.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_15FJforGLM.txt"),"w")
	cmd = "{MACHETE}/parse_FJ_ID_for_GLM.sh {OUTPUT_DIR} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


###### ESTIMATE LIGATION ARTIFACT
## Ligation artifact refers to the rate at which cDNA are artificially ligated, producing what appears to be a true junctional alignment.  In this step, reads that remained unaligned to the Far Junctions Bowtie index that are located in FarJuncSecondary are then aligned to new indel fasta indices under the assumption that in any local area, the rate of false ligation of a junction is similar to the rate of false ligation of any two cDNAs that do not end at exon boundaries.
#
###
### MakeIndelFiles.sh is a shell script that calls the python script AddIndelsToFasta.py.  It takes the FarJunctions fasta files as inputs (FJDir/fasta/<STEM>_FarJunctions.fa) and outputs five files called FJDir/FarJuncIndels/<STEM>/<STEM>_FJ_Indels_1,2,3,4,5.fa where the numbers 1-5 indicate the N's inserted on each side of the breakpoint or deletions on each side of the breakpoint.  For example, the FJ_indels_3 file is the same as the FarJunctions.fa file except every sequence has 3 N's on each side of the breakpoint (total of 6 N's inserted at the breakpoint), or 3 bases are deleted from each exon on each side of the breakpoint (total of 6 deletions at the breakpoint).
####

#j10_id
print("make indel files")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_10FJIndels.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_10FJIndels.txt"),"w")
	cmd = "{MACHETE}/MakeIndelFiles.sh {OUTPUT_DIR} {NumIndels} {MACHETE} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,NumIndels=NumIndels,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

#
#
# make Bowtie Indices for Far Junc Indel files
## The shell script BowtieIndexFJIndels.sh calls bowtie to index the indels_N.fa files that were created in the previous step.  The indels are output to the directory FJDir/BowtieIndels/<STEM>/<STEM>_Indels_N where N is the number of indels in that index.
#
print("index indels")
for i in range(1,NumIndels + 1):
	processes = {}
	for index in range(1,NUM_FILES + 1):
		stdout = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_out_11indexindels.txt".format(i=i,index=index)),"w")
		stderr = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_err_11indexindels.txt".format(i=i,index=index)),"w")
		cmd = "{MACHETE}/BowtieIndexFJIndels.sh {FarJuncIndelsDir} {i} {BowtieIndelsDir} {OUTPUT_DIR} {index}".format(MACHETE=MACHETE,FarJuncIndelsDir=FarJuncIndelsDir,BowtieIndelsDir=BowtieIndelsDir,OUTPUT_DIR=OUTPUT_DIR,i=i,index=index)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
	checkProcesses(processes)


##Align FarJuncSecondary (unaligned to FJ index) to FJ indels

# This section calls the shell BowtieAlignFJIndels.sh to align the fq files FJdir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_1/2.fq to the Bowtie2 indices of the far junction indels created in the previous step.  Aligned indels are stored in FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam where N is the number of indels in the Bowtie index where the reads aligned.  The bowtie parameters include a max of ~4 mismatches / 100 basepairs, max #N is the read length, and prohibits gapped alignments or gaps in the read.

#
##
print("align to indels")
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
for i in range(1,NumIndels + 1):
	processes = {}
	for index in range(1,NUM_FILES + 1):
		stdout = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_out_12alignindels.txt".format(i=i,index=index)),"w")
		stderr = open(os.path.join(LOG_DIR,"NumIndes{i}_{index}_err_12alignindels.txt".format(i=i,index=index)),"w")
		cmd = "{MACHETE}/BowtieAlignFJIndels.sh {OUTPUT_DIR} \"{BOWTIEPARAMETERS}\" {i} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,BOWTIEPARAMETERS=BOWTIEPARAMETERS,i=i,index=index)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
	checkProcesses(processes)


## Calls FindAlignmentArtifact_SLURM.sh which is a shell that calls MakeIndelsHisto.py.  The MakeIndelsHisto.py script reads in the aligned indels from the sam files FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam.  It concatenates all the indels_1,2,3,4,5.sam files into a larger file All_<STEM>_1/2_indels.sam in the same directory as the original sam files. In the All_<STEM>_1/2_indels.sam files, any indels that did not overlap the junction by the user specified # base pairs around the breakpoint are removed. Additionally, since the FarJuncSecondary files were aligned independently to the indels_1-5 bowtie indices, the same read could align to multiple indices.  Therefore for each read, the read with the best alignment score is placed in the All_<STEM>_1/2_indels sam file, and all other alignments to other indices are discarded.
#  Then the python script checks the FJdir/FarJunctionAlignments/<STEM>/ sam files and creates an array for each junction.  The array is of length 2*#indels+1. In the case of 5 indels, the length is 11 and represents the number reads that aligned to the junction with [5Del, 4Del, 3Del, 2Del, 1Del, aligned to junction exactly, 1Ins, 2Ins, 3Ins, 4Ins, 5Ins]
## the junction name and this array are output to FJDir/IndelsHistogram/indels_<STEM>_1/2.txt.  These outputs will be used to generate the Appended reports later.
#j14_id
print("make indels histo")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_13filterIndels.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_13filterIndels.txt"),"w")
	cmd = "{MACHETE}/FindAlignmentArtifact_SLURM.sh {OUTPUT_DIR} {NUMBASESAROUNDJUNC} {NumIndels} {MACHETE} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,NumIndels=NumIndels,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


###
####
##### REG INDELS ################
##
# Align unaligned files to the expanded reg junctions with indels
AlignedIndels = os.path.join(CIRCPIPE_DIR,"orig/RegIndelAlignments")
if not os.path.exists(AlignedIndels):
	os.makedirs(AlignedIndels)

#  To train the GLM, indel alignments are also created for the linear junctions.  The reference index of indels to the linear junctions is static and has already been created and is referenced above as "REG_INDEL_INDICES" on line 44.  The script AlignRegIndels calls bowtie to align reads that were unaligned the the KNIFE indices (in KNIFEdir/orig/unaligned/*.fq) to the REG_INDEL_INDICES, with the parameters of 1) approx 4 mismatches / 100 bases, maximum number N's = readlength, and no gapped alignments or read gaps.
#j16_id
print("Aligning unaligned files to linear junc indels")
BOWTIEPARAMETERS = "--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
for i in range(1,NumIndels + 1):
	processes = {}
	for index in range(1,NUM_FILES + 1):
	 	stdout = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_out_15AlignRegIndels.txt".format(i=i,index=index)),"w")
		stderr = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_err_15AlignRegIndels.txt".format(i=i,index=index)),"w")
		cmd = "{MACHETE}/AlignUnalignedtoRegIndel.sh {CIRCPIPE_DIR} {i} {OUTPUT_DIR} \"{BOWTIEPARAMETERS}\" {REG_INDEL_INDICES} {index}".format(MACHETE=MACHETE,CIRCPIPE_DIR=CIRCPIPE_DIR,i=i,OUTPUT_DIR=OUTPUT_DIR,BOWTIEPARAMETERS=BOWTIEPARAMETERS,REG_INDEL_INDICES=REG_INDEL_INDICES,index=index)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
	checkProcesses(processes)

###
####
######## MAKE REG AND FJ INDELS CLASS OUTPUT FILES ###########
####
###
### reg indels class output file
#Calls RegIndelsClassID.sh shell which calls RegIndels_ClassIDFile.py.
# Inputs are the regular indel alignments located at KNIFEdir/orig/RegIndelAlignments/<STEM>/unaligned_<STEM>_1/2_indel1,2,3,4,5.sam.  These are concatenated into a single file in the same directory called All_<STEM>_1/2_Regindels.sam.  In the concatenation step, like for the far junctions indels, any reads are omitted if they fail to overlie the junction by the user specified overlap, and also if a read aligns to multiple indel indices, the one with the best alignment score is put into the concatenated file.
## Then, the partner reads of all reads from All_<STEM>_1/2_Regindels.sam are identified and labeled as "good" or "bad".  If a read partner is found in genome, it has priority over transcriptome, which has priority over reg, and finally junc.  A far junction R2 cannot be found in another dictionary, as FJ reads are generated from previously unaligned reads.  If the read partner is in genome, it must be located on the same chromosome, on the opposite reference strand from R1.  If a read partner is in reg, then the downstream exon must be upstream of the uptstream reg indel exon, or the upstream read partner exon must be downstream of the downstream reg indel exon, on the same chromosome. Reference strands must be opposite.    If the partner is in junc, then the reg indel alignment must be located inside the circle created by the scrambled junction, and on the opposite reference strand.  In this manner, class input files are generated for the reg indels, which are located at KNIFE dir/circReads/ids/<STEM>_output_RegIndel.txt
#

#j18_id
print("Reg Indels Class Output")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_18RegIndelsClassOutput.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_18RegIndelsClassOutput.txt"),"w")
	cmd = "{MACHETE}/RegIndelsClassID.sh {OUTPUT_DIR} {CIRCPIPE_DIR} {NUMBASESAROUNDJUNC} {MACHETE} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

# FJ indels class output file
## This script calls FJIndelsClassID.sh
## This takes the FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/All_<STEM>_1/2_FJindels.sam and identifies read partners.  The same criteria to identify read partners as FarJuncNaiveReport.sh are applied (see above).
## Output files are placed into FJDir/GLM_classInput/<STEM>_output_FJIndels.txt
#j19_id
print("FJ Indels Class Output")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = open(os.path.join(LOG_DIR,str(index) + "_out_19FJIndelsClassOutput.txt"),"w")
	stderr = open(os.path.join(LOG_DIR,str(index) + "_err_19FJIndelsClassOutput.txt"),"w")
	cmd = "{MACHETE}/FJIndelsClassID.sh {OUTPUT_DIR} {CIRCPIPE_DIR} {NUMBASESAROUNDJUNC} {MACHETE} {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
	print(cmd)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
print("FJ Indels Class Output")


