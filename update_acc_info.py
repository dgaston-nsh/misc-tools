#!/usr/bin/env python

import csv
import sys
import argparse

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--case', help="Current case tracking file")
    parser.add_argument('-o', '--output', help='Output file name for updated case tracking file')
    parser.add_argument('-u', '--update', help="Update file to add")
    parser.add_argument('-t', '--type', help="Type of data update (WGS or WTS)")

    args = parser.parse_args()
    args.logLevel = "INFO"

    case_data = defaultdict(lambda: dict)

    with open(args.case, 'r') as input:
        reader = csv.DictReader(input, delimiter=',')
        for row in reader:
            case_data[row['ACC ID']['External ID'] = row['External ID']
            case_data[row['ACC ID']['Cohort'] = row['Cohort']
            case_data[row['ACC ID']['Status'] = row['Status']
            case_data[row['ACC ID']['RNA Status'] = row['RNA Status']
            case_data[row['ACC ID']['DNA Status'] = row['DNA Status']
            case_data[row['ACC ID']['# Reads RNA (Millions)'] = row['# Reads RNA (Millions)']
            case_data[row['ACC ID']['Normal Coverage'] = row['Normal Coverage']
            case_data[row['ACC ID']['Tumour Coverage'] = row['Tumour Coverage']
            case_data[row['ACC ID']['Sequencing Centre'] = row['Sequencing Centre']
            case_data[row['ACC ID']['Notes'] = row['Notes']

    with open(args.update, 'r') as input:
        reader = csv.DictReader(delimiter=',')
        for row in reader:
            if args.type == 'WGS':
                if row['ACC ID'] in case_data:
                    # Add data
                else:
                    case_data[row['ACC ID']['Cohort'] = row['Cohort']

                    if row['Sample type'] = 'Tumor':
                        case_data[row['ACC ID']['Tumour Coverage'] = float(row['BWA_Aligned_dedup_cov'])
                        if float(row['BWA_Aligned_dedup_cov']) > 80.0:

                    elif row['Sample type'] = 'Normal':
                        case_data[row['ACC ID']['Normal Coverage'] = float(row['BWA_Aligned_dedup_cov'])
                    else:
                        sys.stderr.write(f"ERROR: Unknown sample type ({row['Sample type']})\n")
                        sys.exit()
