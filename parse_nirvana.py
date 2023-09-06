#!/usr/bin/env python

import os
import sys
import gzip
import json
import argparse

import nirvana

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input JSON file to read")
    parser.add_argument('-o', '--output', help="Output filename")
    args = parser.parse_args()
    args.logLevel = "INFO"

    header = "SampleID\tVarID\tGene\tg.HGVS\tc.HGVS\tp.HGVS\t"

    data = nirvana.parseNirvana(args.input)
    with open(args.output, 'w') as outfile:
        outfile.write(f"{header}\n")
        for position in data['positions']:
            if variants_field in position_dict:
                for var_dict in position_dict[variants_field]:
                    if 'clinvar' in var_dict:
                        #clinvar stuff
                    if 'cosmic' in var_dict:
                        #cosmic stuff
                    if 'gnomad' in var_dict:
                        #GnomAD stuff
                    varid = f"{position['chromosome']}-{position['position']}-{position['refAllele']}-{var_dict['altAllele']}"
                    outfile.write(f"{data['samples'][0]}\t")
                    outfile.write(f"{varid}\t")
                    outfile.write(f"{var_dict['']}\t")
                    outfile.write(f"{var_dict['hgvsg']}\t")
