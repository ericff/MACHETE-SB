import re
import os
import subprocess
import csv

# to call this:
# cmd="python {INSTALLDIR}/parse_to_remove_FJ.py --stem {STEM} --outputdir {OUTPUT_DIR}".format(INSTALLDIR=INSTALLDIR, STEM=STEM, OUTPUT_DIR=OUTPUT_DIR)


# python parse_to_remove_FJ.py --stem 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--stem", help="STEM string")
parser.add_argument("--outputdir", help="OUTPUT_DIR")



args = parser.parse_args()

# first get STEM, assume it is passed from run.py
STEM = args.stem
outputdir = args.outputdir


naive_report = os.path.join(outputdir,"reports",STEM + "_naive_report.txt")

fastadir = os.path.join(outputdir,"fasta")

rawcsv = csv.reader(open(naive_report), delimiter='\t')
first_column= []
for row in rawcsv:
    first_column.append(row[0])
# drop first element as it's from the header:
first_column.pop(0)

first_col_with_arrow = [(">" + x) for x in first_column]


fasta_stem_dir = os.path.join(fastadir,STEM)

posschrfjfiles = os.listdir(fasta_stem_dir)

chrfjfiles = [os.path.join(fasta_stem_dir,x) for x in posschrfjfiles if (re.search(pattern='chr.*FarJunctions.fa', string=x))]


newfile = os.path.join(fasta_stem_dir,"temp211.txt")
for thisfile in chrfjfiles:
    # clear newfile in case already used
    open(newfile, 'w').close()
    
    with open(thisfile, 'r') as ff:
        lines = ff.read().splitlines()
    for index, thisline in enumerate(lines):
        if (thisline in first_col_with_arrow):
            with open(newfile, 'a') as nn:
                nn.write(thisline +"\n")
                if (index < len(lines)):
                    nn.write(lines[(index+1)] +"\n")
    # then write new file over the current file
    cmd = "mv -f " + newfile + " " + thisfile
    print(cmd)
    subprocess.call(cmd, shell=True)
    
