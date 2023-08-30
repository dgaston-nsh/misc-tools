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
            elements = d.split('_')
            new_dir = os.path.join(rootdir, elements[1])
            sys.stdout.write(f"Renaming Directory {d} to {new_dir}\n")
            os.rename(d, new_dir)
            os.chdir(new_dir)
            fastqs = glob.glob('*_R1_001.fastq.gz')
            for fastq1 in fastqs:
                elements = fastq1.split('_')
                sample_id = elements[1]

                fastq2 = f"{elements[0]}_{elements[1]}_{elements[2]}_{elements[3]}_{elements[4]}_R2_001.fastq.gz"
                # index1 = f"{elements3[0]}_I1.fastq.gz"
                # index2 = f"{elements3[0]}_I2.fastq.gz"

                sys.stdout.write(f"Renaming {fastq1} to {sample_id}.R1.fastq.gz\n")
                os.rename(fastq1, f"{sample_id}.R1.fastq.gz")

                sys.stdout.write(f"Renaming {fastq2} to {sample_id}.R2.fastq.gz\n")
                os.rename(fastq2, f"{sample_id}.R2.fastq.gz")

                # sys.stdout.write(f"Renaming {index1} to {sample_id}.I1.fastq.gz\n")
                # os.rename(index1, f"{sample_id}.I1.fastq.gz")

                # sys.stdout.write(f"Renaming {index2} to {sample_id}.I2.fastq.gz\n")
                # os.rename(index2, f"{sample_id}.I2.fastq.gz")
            os.chdir(rootdir)
