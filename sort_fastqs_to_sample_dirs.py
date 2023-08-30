#!/usr/bin/env python

# Standard packages
import os
import glob
import shutil

if __name__ == "__main__":
    cwd = os.getcwd()
    fastqs = glob.glob('*_R1.fastq.gz')
    for fastq1 in fastqs:
        elements = fastq1.split('_R1')
        sample_id = elements[0]
        fastq2 = f"{sample_id}_R2.fastq.gz"
        # index1 = f"{sample_id}_I1.fastq.gz"
        # index2 = f"{sample_id}_I2.fastq.gz"

        os.mkdir(sample_id)
        shutil.move(fastq1, f"{sample_id}/{fastq1}")
        shutil.move(fastq2, f"{sample_id}/{fastq2}")
        # shutil.move(index1, f"{sample_id}/{index1}")
        # shutil.move(index2, f"{sample_id}/{index2}")
