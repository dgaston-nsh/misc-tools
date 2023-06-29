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
            elements = d.split('.')
            sys.stdout.write(f"Renaming Directory {d} to {elements[-1]}")
            os.chdir(d)
            fastqs = glob.glob('*_R1.fastq.gz')
            for fastq1 in fastqs:
                elements = fastq1.split('.')
                fastq_id = elements[-3]
                elements2 = fastq_id.split('_')
                sample_id = elements2[0]

                elements3 = fastq1.split('_R1')
                fastq2 = f"{elements3[0]}_R2.fastq.gz"
                index1 = f"{elements3[0]}_I1.fastq.gz"
                index2 = f"{elements3[0]}_I2.fastq.gz"

                #sys.stdout.write(f"Renaming {fastq1} to {sample_id}.R1.fastq.gz\n")
                #os.rename(fastq1, f"{sample_id}.R1.fastq.gz")

                #sys.stdout.write(f"Renaming {fastq2} to {sample_id}.R2.fastq.gz\n")
                #os.rename(fastq2, f"{sample_id}.R2.fastq.gz")

                #sys.stdout.write(f"Renaming {index1} to {sample_id}.I1.fastq.gz\n")
                #os.rename(index1, f"{sample_id}.I1.fastq.gz")

                #sys.stdout.write(f"Renaming {index2} to {sample_id}.I2.fastq.gz\n")
                #os.rename(index2, f"{sample_id}.I2.fastq.gz")
            os.chdir(rootdir)
