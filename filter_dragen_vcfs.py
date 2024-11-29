#!/usr/bin/env python
import os
import sys
import argparse

from ctvcf2 import VCF

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', 'Input list of samples')

    args = parser.parse_args()
    args.logLevel = "INFO"

    with open(args.samples, 'r') as input:
        reader = csv.DictReader(input, delimiter=',', quotechar='"')
        for row in reader:
            sample_output_dir = os.path.join(os.getcwd(), f"{row['RGSM']}_wgs")
            vcf = os.path.join(sample_output_dir, f"{row['RGSM']}_wgs.hard-filtered.vcf.gz")

            with open(f"{row['RGSM']}_wgs.hard-filtered.csv", w) as outfile:
                for variant in VCF(vcf):
