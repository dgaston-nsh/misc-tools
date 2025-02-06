#!/usr/bin/env python

# Import packages
import os
import csv
import sys
import gzip
import json
import argparse
import resource
import xlsxwriter

from collections import defaultdict

def get_json_variants_fname(sample_id, sample_dir):
    annotated_json = os.path.join(sample_dir, f"{sample_id}.hard-filtered.vcf.annotated.json.gz")

    return annotated_json

def get_output_csv_file(sample_id, output_dir):
    filename = os.path.join(output_dir, f"{sample_id}.filtered_variants.csv")

    return filename

def get_output_header():
    header = ['VarID', 'Chr', 'Pos', 'Ref', 'Alt', 'Gene', 'TranscriptID', 'c.HGVS', 'p.HGVS',
    'biotype', 'ProteinID', 'CovDepth', 'Filters', 'RefDepth', 'AltDepth', 'dbSNP',
    'COSMIC_IDs', 'ClinVar_IDs', 'Clingen_IDs',
    'gnomAD_allAF', 'gnomAD_maleAf', 'gnomAD_femaleAf', 'gnomAD_afrAf', 'gnomAD_amrAf',
    'gnomAD_easAf', 'gnomAD_finAf', 'gnomAD_nfeAf', 'gnomAD_asjAf', 'gnomAD_othAf',
    'gnomAD_allAn', 'gnomAD_maleAn', 'gnomAD_femaleAn', 'gnomAD_afrAn', 'gnomAD_amrAn',
    'gnomAD_easAn', 'gnomAD_finAn', 'gnomAD_nfeAn', 'gnomAD_asjAn', 'gnomAD_othAn',
    'gnomAD_allAc','gnomAD_maleAc', 'gnomAD_femaleAc', 'gnomAD_afrAc', 'gnomAD_amrAc',
    'gnomAD_easAc', 'gnomAD_finAc', 'gnomAD_nfeAc', 'gnomAD_asjAc', 'gnomAD_othAc',
    'gnomAD_allHc', 'gnomAD_afrHc', 'gnomAD_amrHc', 'gnomAD_easHc', 'gnomAD_finHc',
    'gnomAD_nfeHc', 'gnomAD_asjHc', 'gnomAD_othHc', 'gnomAD_maleHc', 'gnomAD_femaleHc',
    'gnomAD_controlsAllAf', 'gnomAD_controlsAllAn', 'gnomAD_controlsAllAc',
    'gnomAD_lowComplexityRegion', 'gnomAD_failedFilter',
    'topmed_allAc', 'topmed_allAn', 'topmed_allAf', 'topmed_allHc', 'topmed_failedFilter',
    'PrimateAI_ScorePercentile', 'PhyloP', 'PolyPhenScore', 'PolyPhenPred', 'SiftScore',
    'SiftPred', 'REVEL', 'SpliceAI_accGainScore', 'SpliceAI_accLossScore',
    'SpliceAI_donGainScore', 'SpliceAI_donLossScore', 'g.HGVS']

    return header

def parseNirvana(sample_id, file, output_fname):
    output_header = get_output_header()

    is_header_line = True
    is_position_line = False
    is_gene_line = False

    gene_section_line = '],"genes":['
    end_line = ']}'

    data = dict()

    data['header'] = ''
    data['samples'] = []
    data['positions'] = []
    data['genes'] = []

    with gzip.open(file, 'rt') as f:
        data['position_count'] = 0
        data['gene_count'] = 0

        for line in f:
            trimmed_line = line.strip()
            if is_header_line:
                data['header'] = trimmed_line[10:-14]
                header_dict = json.loads(data['header'])
                data['samples'] = header_dict['samples']
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
                    position_dict = json.loads(trimmed_line.rstrip(','))
                    # data['positions'].append(position_dict)
                    data['position_count'] += 1

                    parse_nirvana_position_dict(position_dict)
                if is_gene_line:
                    data['genes'].append(trimmed_line.rstrip(','))
                    data['gene_count'] += 1

def parse_nirvana_position_dict(position_dict):
    output_lines = list()

    return output_lines

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples', help="Input file for samples")

    args = parser.parse_args()
    args.logLevel = "INFO"

    main_dir = os.getcwd()
    output_dir = os.path.join(main_dir, "ParsedVariantReports")

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    with open(args.samples, 'r') as input:
        reader = csv.DictReader(input, delimiter=',', quotechar='"')
        for row in reader:
            sample_output_dir = os.path.join(main_dir, f"{row['sampleID']}_wgs")

            json_fname = get_json_variants_fname(row['sampleID'], sample_output_dir)
            output_fname = get_output_csv_file(row['sampleID'], sample_output_dir)
            sys.stdout.write(f"Parsing {json_fname}\n")

            parseNirvana(sample_id, json_fname, output_fname)
