#!/usr/bin/env python

import csv
import sys
import argparse

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--case_data', help="Input file name")
    parser.add_argument('-e', '--acc_and_external_ids', help="Output file name")
    parser.add_argument('-o', '--output', help='output file name')

    args = parser.parse_args()
    args.logLevel = "INFO"

    acc_external_map = dict()

    with open(args.acc_and_external_ids, 'r') as external:
        reader = csv.reader(external, delimiter = ',')
        next(reader)
        for row in reader:
            acc_external_map[row[13]] = row[11]

    with open(args.case_data, 'r') as input:
        with open(args.output, 'w') as output:
            reader = csv.reader(input, delimiter=',')
            writer = csv.writer(output, delimiter=',')
            writer.writerow(['ACC ID', 'External ID', 'Cohort', 'Status', 'RNA Status', 'DNA Status', '# Reads RNA (Millions)', 'Normal Coverage', 'Tumour Coverage', 'Sequencing Centre', 'Notes'])
            next(reader)

            for row in reader:
                try:
                    external_id = acc_external_map[row[0]]
                except:
                    external_id = ''

                writer.writerow([row[0], external_id, row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10]])
