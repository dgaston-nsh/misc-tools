#!/usr/bin/env python

import os
import sys
import gzip
import json
import argparse

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input JSON file to read")
    args = parser.parse_args()
    args.logLevel = "INFO"

    header = ''
    positions = []
    genes = []

    is_header_line = True
    is_position_line = False
    is_gene_line = False

    gene_section_line = '],"genes":['
    end_line = ']}'

    with gzip.open(args.input, 'rt') as f:
        position_count = 0
        gene_count = 0

        for line in f:
            trimmed_line = line.strip()
            if is_header_line:
                header = trimmed_line[10:-14]
                is_header_line = False
                is_position_line = True
                continue
            if trimmed_line == gene_section_line:
                is_gene_line = True
                is_position_line = False
                continue
            elif trimmed_line == end_line:
                break
            else:
                if is_position_line:
                    positions.append(trimmed_line.rstrip(','))
                    position_count += 1
                if is_gene_line:
                    genes.append(trimmed_line.rstrip(','))
                    gene_count += 1

    print ('header object:', header)
    print ('number of positions:', position_count)
    print ('number of genes:', gene_count)
