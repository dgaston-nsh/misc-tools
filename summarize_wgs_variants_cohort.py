#!/usr/bin/env python

# Import packages
import os
import csv
import sys
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

    main_dir = os.getcwd()
