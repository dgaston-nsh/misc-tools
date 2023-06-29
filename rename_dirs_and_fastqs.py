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
            fastqs = glob.glob('*_R1.fastq.gz')
            for fastq1 in fastqs:
                elements = fastq1.split('.')
                fastq_id = elements[-3]
                elements2 = fastq_id.split('_')
                sample_id = elements2[0]
                sys.stdout.write(f"Renaming {fastq1} to {sample_id}.R1.fastq.gz\n")
