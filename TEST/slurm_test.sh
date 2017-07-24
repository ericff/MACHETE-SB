#!/bin/bash -eu

module load machete/current python/2.7.5 R/3.2.2 bowtie2/2.2.8


python /scratch/users/nathankw/MACHETE/run.py /scratch/users/glhsieh/fetalcircpipe/ \
 /scratch/users/nathankw/machete_runs/run1 \
 100000 \
 HG19 \
 8 \
 5
