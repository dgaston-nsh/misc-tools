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
    'biotype', 'Consequence', 'ProteinID', 'CovDepth', 'Filters', 'RefDepth', 'AltDepth', 'VAF', 'dbSNP',
    'COSMIC_IDs', 'COSMIC_NumSamples', 'ClinVar_IDs', 'ClinVarSig', 'Clingen_IDs',
    'gnomAD_allAF', 'gnomAD_maleAf', 'gnomAD_femaleAf', 'gnomAD_afrAf', 'gnomAD_amrAf',
    'gnomAD_easAf', 'gnomAD_finAf', 'gnomAD_nfeAf', 'gnomAD_asjAf', 'gnomAD_othAf',
    'gnomAD_controlsAllAf', 'gnomAD_largestAF', 'topmed_allAf', 'PrimateAI_ScorePercentile', 'PhyloP',
    'PolyPhenScore', 'PolyPhenPred', 'SiftScore', 'SiftPred', 'REVEL',
    'SpliceAI_accGainScore', 'SpliceAI_accLossScore', 'SpliceAI_donGainScore', 'SpliceAI_donLossScore', 'g.HGVS']

    return header

def parseBasicInfo(out, position_dict,samples_dict):
    out['Chr'] = position_dict['chromosome']
    out['Pos'] = position_dict['position']
    out['Ref'] = position_dict['refAllele']
    out['Filters'] = position_dict['filters'][0]

    if 'totalDepth' in samples_dict:
        out['CovDepth'] = samples_dict['totalDepth']

    if 'alleleDepths' in samples_dict:
                    out['RefDepth'] = samples_dict['alleleDepths'][0]

    return out

def parseBasicVarInfo(out, var_dict, samples_dict, var_index, allele_index):
    out['Alt'] = var_dict['altAllele']
    out['VarID'] = var_dict['vid']

    if 'hgvsg' in var_dict:
        out['g.HGVS'] = var_dict['hgvsg']

    if 'variantFrequencies' in samples_dict:
        out['VAF'] = float(samples_dict['variantFrequencies'][var_index])

    if 'alleleDepths' in samples_dict:
        out['AltDepth'] = int(samples_dict['alleleDepths'][allele_index])

    if 'phylopScore' in var_dict:
        out['PhyloP'] = var_dict['phylopScore']

    if 'dbsnp' in var_dict:
        out['dbSNP'] = ','.join(var_dict['dbsnp'])

    if 'PrimateAI-3D' in var_dict:
        primateai_dict = var_dict['PrimateAI-3D']
        out['PrimateAI_ScorePercentile'] = primateai_dict['scorePercentile']

    if 'revel' in var_dict:
        revel_dict = var_dict['revel']
        out['REVEL'] = revel_dict['score']

    out = parseClinVar(var_dict, out)
    out = parseCOSMIC(var_dict, out)
    out = parsePopAFs(var_dict, out)
    out = parseSpliceAI(var_dict, out)

    return out

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

    return out

def parsePopAFs(var_dict, out):
    gnomad_populations = ["allAf", "afrAf", "amrAf", "easAf", "finAf", "nfeAf", "asjAf", "sasAf", "othAf", "controlsAllAf"]

    if 'gnomad' in var_dict:
        gnomad_dict = var_dict['gnomad']
        out['gnomAD_allAF'] = gnomad_dict['allAf']
        largest_af = 0.0
        for pop in gnomad_populations:
            try:
                freq = float(gnomad_dict.get(pop))
            except:
                freq = 0
            out[f"gnomAD_{pop}"] = freq

            if freq > largest_af:
                largest_af = freq
        out['gnomAD_largestAF'] = largest_af

    if 'topmed' in var_dict:
        topmed_dict = var_dict['topmed']
        out['TopMed_allAF'] = float(topmed_dict['allAf'])

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

def parseTranscriptInfo(out,transcript_dict):
    out['Gene'] = transcript_dict['hgnc']
    out['biotype'] = transcript_dict['bioType']
    out['Consequence'] = ','.join(transcript_dict['consequence'])
    out['TranscriptID'] = transcript_dict['transcript']

    if 'proteinId' in transcript_dict:
        out['ProteinID'] = transcript_dict['proteinId']
    else:
        out['ProteinID'] = "-"

    if 'hgvsc' in transcript_dict:
        out['c.HGVS'] = transcript_dict['hgvsc']
    else:
        out['c.HGVS'] = "-"

    if 'hgvsp' in transcript_dict:
        out['p.HGVS'] = transcript_dict['hgvsp']
    else:
        out['p.HGVS'] = "-"

    if 'polyPhenScore' in transcript_dict:
        out['polyPhenScore'] = transcript_dict['polyPhenScore']
    else:
        out['polyPhenScore'] = '-'

    if 'polyPhenPred' in transcript_dict:
        out['polyPhenPred'] = transcript_dict['polyPhenPred']
    else:
        out['polyPhenPred'] = '-'

    if 'siftScore' in transcript_dict:
        out['siftScore'] = transcript_dict['siftScore']
    else:
        out['siftScore'] = '-'

    if 'siftPred' in transcript_dict:
        out['siftPred'] = transcript_dict['siftPrediction']
    else:
        out['siftPred'] = '-'

    return out

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
    # data['positions'] = []
    data['genes'] = []

    with open(output_fname, 'w') as output:
        writer = csv.DictWriter(output, fieldnames=output_header)
        writer.writeheader()
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

                        if len(position_dict['filters']) == 1:
                            if position_dict['filters'][0] == 'PASS':
                                out = defaultdict(None)
                                var_index = 0
                                allele_index = 1
                                samples_dict = position_dict['samples'][0]

                                out = parseBasicInfo(out, position_dict,samples_dict)

                                # Set a minimum coverage of 5X
                                if out['CovDepth'] >= 5:
                                    if 'variants' in position_dict:
                                        for var_dict in position_dict['variants']:
                                            out = parseBasicVarInfo(out, var_dict, samples_dict, var_index, allele_index)

                                            # Set VAF of 0.35 threshold
                                            if out.get('VAF') >= 0.35:
                                                # Set the 1% Allele Frequency cutoff here using the Controls
                                                if out.get('gnomAD_controlsAF') >= 0.005:
                                                    if 'transcripts' in var_dict:
                                                        for transcript_dict in var_dict['transcripts']:
                                                            out = parseTranscriptInfo(out, transcript_dict)
                                                            writer.write(out)
                                                    else:
                                                        writer.write(out)

                                            var_index += 1
                                            allele_index += 1
                                else:
                                    sys.stderr.write(f"Position in Annotated JSON with no Variant entries: {out['Chr']}-{out['Pos']}-{out['Ref']}\n")



                    if is_gene_line:
                        data['genes'].append(trimmed_line.rstrip(','))
                        data['gene_count'] += 1

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

            parseNirvana(row['sampleID'], json_fname, output_fname)
