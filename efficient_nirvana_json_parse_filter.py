#!/usr/bin/env python

import csv
import sys
import gzip
import json

import numpy as np

from collections import defaultdict

def parseClinVar(var_dict, out):
    if 'clinvar' in var_dict:
        clinvar_ids = []
        clinvar_sigs = []
        exact_match = False
        for clinvar_dict in var_dict['clinvar']:
            clinvar_ids.append(clinvar_dict['id'])
            clinvar_sigs.append("|".join(clinvar_dict['significance']))
            if clinvar_dict['refAllele'] == out['Ref']:
                if clinvar_dict['altAllele'] == out['Alt']:
                    exact_match = True
                    out['ClinVarID'] = clinvar_dict['id']
                    out['ClinVarSig'] = "|".join(clinvar_dict['significance'])
        if not exact_match:
            ids_string = ",".join(clinvar_ids)
            sigs_string = ",".join(clinvar_sigs)
            out['ClinVarID'] = f"Overlapping ({ids_string})"
            out['ClinVarSig'] = f"Overlapping ({sigs_string})"
    else:
        out['ClinVarID'] = "-"
        out['ClinVarSig'] = "-"

    return out

def parseCOSMIC(var_dict, out):
    if 'cosmic' in var_dict:
        cosmic_ids = []
        total_num_overlap = 0
        exact_match = False
        for cosmic_dict in var_dict['cosmic']:
            cosmic_ids.append(cosmic_dict['id'])
            total_num_overlap += int(cosmic_dict['numSamples'])
            if cosmic_dict['refAllele'] == out['Ref']:
                if cosmic_dict['altAllele'] == out['Alt']:
                    exact_match = True
                    out['COSMIC_ID'] = cosmic_dict['id']
                    out['COSMIC_NumSamples'] = int(cosmic_dict['numSamples'])
        if not exact_match:
            ids_string = ",".join(cosmic_ids)
            out['COSMIC_ID'] = f"Overlapping ({ids_string})"
            out['COSMIC_NumSamples'] = total_num_overlap
    else:
        out['COSMIC_ID'] = "-"
        out['COSMIC_NumSamples'] = 0

    return out

def parseSpliceAI(var_dict, out):
    if 'spliceAI' in var_dict:
        spliceai_list = var_dict['spliceAI']
        spliceai_dict = spliceai_list[0]

        if 'acceptorGainScore' in spliceai_dict:
            out['SpliceAI_accGainScore'] = float(spliceai_dict['acceptorGainScore'])
        else:
            out['SpliceAI_accGainScore'] = 0.0

        if 'acceptorLossScore' in spliceai_dict:
            out['SpliceAI_accLossScore'] = float(spliceai_dict['acceptorLossScore'])
        else:
            out['SpliceAI_accLossScore'] = 0.0

        if 'donorGainScore' in spliceai_dict:
            out['SpliceAI_donGainScore'] = float(spliceai_dict['donorGainScore'])
        else:
            out['SpliceAI_donGainScore'] = 0.0

        if 'donorLossScore' in spliceai_dict:
            out['SpliceAI_donLossScore'] = float(spliceai_dict['donorLossScore'])
        else:
            out['SpliceAI_donLossScore'] = 0.0
    else:
        out['SpliceAI_accGainScore'] = 0.0
        out['SpliceAI_accLossScore'] = 0.0
        out['SpliceAI_donGainScore'] = 0.0
        out['SpliceAI_donLossScore'] = 0.0

    return out

def parsePopAFs(var_dict):
    afs = dict()
    gnomad_populations = ["allAf", "afrAf", "amrAf", "easAf", "finAf", "nfeAf", "asjAf", "sasAf", "othAf"]

    if 'gnomad' in var_dict:
        gnomad_dict = var_dict['gnomad']
        afs['gnomAD_allAF'] = gnomad_dict['allAf']
        largest_af = 0.0
        for pop in gnomad_populations:
            try:
                freq = float(gnomad_dict.get(pop))
            except:
                freq = 0
            if freq > largest_af:
                largest_af = freq
        afs['gnomAD_largestAF'] = largest_af
    else:
        afs['gnomAD_allAF'] = 0.0
        afs['gnomAD_largestAF'] = 0.0

    if 'topmed' in var_dict:
        topmed_dict = var_dict['topmed']
        afs['TopMed_allAF'] = float(topmed_dict['allAf'])
    else:
        afs['TopMed_allAF'] = 0.0

    return afs

def classify_filter(position_dict):
    output = dict()

    var_index = 0
    allele_index = 1

    samples_dict = position_dict['samples'][0]
    depth = int(samples_dict.get('totalDepth'))

    if depth > 5:
        if 'variants' in position_dict:
            for var_dict in position_dict['variants']:
                if 'variantFrequencies' in samples_dict:
                    vaf = float(samples_dict['variantFrequencies'][var_index])
                    if vaf >= 0.35:
                        pop_af_dict = parsePopAFs(var_dict)
                        if pop_af_dict['gnomAD_allAF'] < 0.005 and pop_af_dict['TopMed_allAF'] < 0.005:
                            output['VarID'] = var_dict['vid']
                            output['Depth'] = depth
                            output['VAF'] = vaf
                            output['Filter'] = "|".join(position_dict['filters'])
                            output['gnomAD_allAF'] = pop_af_dict['gnomAD_allAF']
                            output['TopMed_allAF'] = pop_af_dict['TopMed_allAF']
                            output['g.HGVS'] = var_dict['hgvsg']

                            if 'phylopScore' in var_dict:
                                output['PhyloP'] = var_dict['phylopScore']
                            else:
                                output['PhyloP'] = "-"

                            if 'dbsnp' in var_dict:
                                output['dbSNP'] = ','.join(var_dict['dbsnp'])
                            else:
                                output['dbSNP'] = "-"

                            if 'PrimateAI-3D' in var_dict:
                                primateai_dict = var_dict['PrimateAI-3D']
                                output['PrimateAI_ScorePercentile'] = primateai_dict['scorePercentile']
                            else:
                                output['PrimateAI_ScorePercentile'] = '-'

                            if 'revel' in var_dict:
                                revel_dict = var_dict['revel']
                                output['REVEL'] = revel_dict['score']
                            else:
                                output['REVEL'] = '-'

                            output = parseClinVar(var_dict, out)
                            output = parseCOSMIC(var_dict, out)
                            output = parseSpliceAI(var_dict, out)

                            if 'transcripts' in var_dict:
                                for transcript_dict in var_dict['transcripts']:

                var_index += 1
                allele_index += 1

    return output

def parse_json(json_file, output_csv):
    is_header_line = True
    is_position_line = False
    is_gene_line = False

    header = ''
    genes = []

    gene_section_line = '],"genes":['
    end_line = ']}'

    fieldnames = []

    with open(output_csv, w) as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writerheader()
        with gzip.open(file, 'rt') as f:
        position_count = 0
        gene_count = 0

        for line in f:
            trimmed_line = line.strip()
            if is_header_line:
                header = trimmed_line[10:-14]
                header_dict = json.loads(header)
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
                    output = classify_filter(position_dict)
                    if len(output) > 0:
                        writer.writerow(output)
                if is_gene_line:
                    genes.append(trimmed_line.rstrip(','))
                    gene_count += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples', help="Input list of samples")

    args = parser.parse_args()
    args.logLevel = "INFO"

    with open(args.samples, 'r') as input:
        reader = csv.DictReader(input, delimiter=',', quotechar='"')
        for row in reader:
            sample_output_dir = os.path.join(os.getcwd(), f"{row['RGSM']}_wgs")
            json_file = os.path.join(sample_output_dir, f"{row['RGSM']}.hard-filtered.vcf.annotated.json.gz")
            output_csv = os.path.join(sample_output_dir, f"{row['RGSM']}_wgs.hard-filtered.csv"")

            parse_json(json_file, output_csv)
