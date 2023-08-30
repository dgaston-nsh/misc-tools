#!/usr/bin/env python

# Standard packages
import os
import sys
import glob
import shutil

if __name__ == "__main__":
    rootdir = os.getcwd()
    for file in os.listdir(rootdir):
        d = os.path.join(rootdir, file)
        if os.path.isdir(d):
            os.chdir(d)
            fastqs = glob.glob('*.R1.fastq.gz')
            for fastq1 in fastqs:
                elements = fastq1.split('.')
                sample_id = elements[0]
                fastq2 = f"{sample_id}.R2.fastq.gz"

                sys.stdout.write(f"Renaming {fastq1} to {sample_id}.R1.fastq.gz\n")
                os.rename(fastq1, f"{sample_id}-DNA.R1.fastq.gz")

                sys.stdout.write(f"Renaming {fastq2} to {sample_id}.R2.fastq.gz\n")
                os.rename(fastq2, f"{sample_id}-DNA.R2.fastq.gz")
