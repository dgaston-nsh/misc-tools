#!/usr/bin/env python

# Import packages
import os
import csv
import sys
import glob
import argparse

from collections import defaultdict

# header = ['VarID', 'Chr', 'Pos', 'Ref', 'Alt', 'genotype', 'genotypeQuality',
# 'Gene and Transcript Info', 'CovDepth', 'Filters', 'RefDepth', 'AltDepth', 'VAF', 'dbSNP', 'COSMIC_ID',
# 'COSMIC_NumSamples', 'ClinVarID', 'ClinVarSig', 'Clingen_IDs',
# 'gnomAD_allAf', 'gnomAD_maleAf', 'gnomAD_femaleAf', 'gnomAD_afrAf', 'gnomAD_amrAf',
# 'gnomAD_easAf', 'gnomAD_sasAf', 'gnomAD_finAf', 'gnomAD_nfeAf', 'gnomAD_asjAf', 'gnomAD_othAf',
# 'gnomAD_controlsAllAf', 'gnomAD_largestAF', 'TopMed_allAF', 'PrimateAI_ScorePercentile', 'PhyloP',
# 'polyPhenScore', 'polyPhenPred', 'siftScore', 'siftPred', 'REVEL',
# 'SpliceAI_accGainScore', 'SpliceAI_accLossScore', 'SpliceAI_donGainScore', 'SpliceAI_donLossScore', 'g.HGVS']

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help="Output file name for variants summary")

    args = parser.parse_args()
    args.logLevel = "INFO"

    variant_counts = defaultdict(list)
    csv_files = glob.glob('*.csv')

    for file in csv_files:
        with open(file, 'r') as csv_file:
            elements = csv_file.split('.')
            sample_id = elements[0]
            reader = csv.reader(csv_file)
            for row in reader:
                variant_counts[row[0]].append(sample_id)

    with open(args.output, 'w') as output:
        output.write("VarID\tNum Instances\tSamples\n")
        for varid in variant_counts:
            count = len(variant_counts[varid])
            samples_string = "|".join(variant_counts[varid])
            output.write(f"{varid}\t{count}\t{samples_string}\n")
