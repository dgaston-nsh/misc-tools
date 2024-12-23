#!/usr/bin/env python

import csv
import sys
import argparse

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input file name")
    parser.add_argument('-o', '--output', help="Output file name")

    args = parser.parse_args()
    args.logLevel = "INFO"

    with open(args.input, 'r') as input:
        with open(args.output, 'w') as output:
            reader = csv.DictReader(input, delimiter=',')
            writer = csv.writer(output, delimiter=',')
            writer.writerow(['ACC ID', 'External ID', 'Cohort', 'Status', 'RNA Status', 'DNA Status', '# Reads RNA (Millions)', 'Normal Coverage', 'Tumour Coverage', 'Sequencing Centre', 'Notes'])
            for row in reader:
                try:
                    rna_reads = int(row['# Million RNA Read Pairs'])
                except:
                    rna_reads = False

                try:
                    norm_cov = int(row['Normal WGS Coverage (30X)'])
                except:
                    norm_cov = False

                try:
                    tumour_cov = int(row['Tumour WGS Coverage (80X/30X)'])
                except:
                    tumour_cov = False

                if norm_cov and tumour_cov:
                    if norm_cov >= 30 and tumour_cov >= 80:
                        dna_status = 'Tier A'
                    elif norm_cov >= 30 and tumour_cov >= 30:
                        dna_status = 'Tier B'
                    elif tumour_cov >= 80:
                        dna_status = 'Incomplete Tier A (Needs Normal Top-Up)'
                    elif tumour_cov >= 30:
                        dna_status = 'Incomplete Tier B (Needs Normal Top-Up)'
                    else:
                        dna_status = 'FAIL'
                elif tumour_cov:
                    dna_status = 'Incomplete - Waiting on Normal'
                elif norm_cov:
                    dna_status = 'Incomplete - Waiting on Tumour'
                else:
                    dna_status = 'Incomplete - No DNA'

                if rna_reads:
                    if rna_reads > 80:
                        rna_status = 'Tier A'
                    else:
                        rna_status = 'Tier B'
                else:
                    rna_status = 'Incomplete - No RNA'

                if dna_status == 'Tier A' and rna_status == 'Tier A':
                    status = 'Tier A'
                elif dna_status == 'Tier A' and rna_status == 'Tier B':
                    status = 'Tier B'
                elif dna_status == 'Tier B' and rna_status == 'Tier A':
                    status = 'Tier B'
                elif rna_status == 'Tier A':
                    if dna_status == 'Incomplete Tier A (Needs Normal Top-Up)':
                        status = 'Incomplete - Tier A'
                    elif dna_status == 'Incomplete Tier B (Needs Normal Top-Up)':
                        status = 'Incomplete - Tier A'
                    else:
                        status = 'Incomplete'
                elif rna_status == 'Tier B':
                    if dna_status == 'Incomplete Tier A (Needs Normal Top-Up)':
                        status = 'Incomplete - Tier B'
                    elif dna_status == 'Incomplete Tier B (Needs Normal Top-Up)':
                        status = 'Incomplete - Tier B'
                    else:
                        status = 'Incomplete'
                elif rna_status == 'Incomplete - No RNA':
                    if dna_status == 'Tier A' or dna_status == 'Tier B':
                        status = "Tier A or B"
                    elif dna_status == 'Incomplete Tier A (Needs Normal Top-Up)':
                        status = 'Incomplete - Tier A or B'
                    elif dna_status == 'Incomplete Tier B (Needs Normal Top-Up)':
                        status = 'Incomplete - Tier B'
                    else:
                        status = ''
                else:
                    status = ''

                writer.writerow([row['ACC-ID'], row['ExternalID'], row['Cohort'], status, rna_status, dna_status, row['# Million RNA Read Pairs'], row['Normal WGS Coverage (30X)'], row['Tumour WGS Coverage (80X/30X)'], row['Sequencing Centre'], ''])
