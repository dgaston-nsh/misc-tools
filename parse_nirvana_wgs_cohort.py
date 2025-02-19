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

from multiprocessing import Pool
from collections import defaultdict

def get_json_variants_fname(sample_id, sample_dir):
    annotated_json = os.path.join(sample_dir, f"{sample_id}_wgs.hard-filtered.vcf.annotated.json.gz")

    return annotated_json

def get_output_csv_file(sample_id, output_dir):
    filename = os.path.join(output_dir, f"{sample_id}.filtered_variants.csv")

    return filename

def get_output_header():
    header = ['VarID', 'Chr', 'Pos', 'Ref', 'Alt', 'genotype', 'genotypeQuality',
    'Gene and Transcript Info', 'CovDepth', 'Filters', 'RefDepth', 'AltDepth', 'VAF', 'dbSNP', 'COSMIC_ID',
    'COSMIC_NumSamples', 'ClinVarID', 'ClinVarSig', 'Clingen_IDs',
    'gnomAD_allAf', 'gnomAD_maleAf', 'gnomAD_femaleAf', 'gnomAD_afrAf', 'gnomAD_amrAf',
    'gnomAD_easAf', 'gnomAD_sasAf', 'gnomAD_finAf', 'gnomAD_nfeAf', 'gnomAD_asjAf', 'gnomAD_othAf',
    'gnomAD_controlsAllAf', 'gnomAD_largestAF', 'TopMed_allAF', 'PrimateAI_ScorePercentile', 'PhyloP',
    'polyPhenScore', 'polyPhenPred', 'siftScore', 'siftPred', 'REVEL',
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

    if 'genotype' in samples_dict:
        out['genotype'] = samples_dict['genotype']

    if 'genotypeQuality' in samples_dict:
        out['genotypeQuality'] = samples_dict['genotypeQuality']

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
        out['gnomAD_allAf'] = gnomad_dict['allAf']
        largest_af = 0.0
        for pop in gnomad_populations:
            try:
                freq = float(gnomad_dict.get(pop))
            except:
                freq = 0.0
            out[f"gnomAD_{pop}"] = freq

            if freq > largest_af:
                largest_af = freq
        out['gnomAD_largestAF'] = largest_af
    else:
        for pop in gnomad_populations:
            out[f"gnomAD_{pop}"] = 0.0

    if 'topmed' in var_dict:
        topmed_dict = var_dict['topmed']
        out['TopMed_allAF'] = float(topmed_dict['allAf'])
    else:
        out['TopMed_allAF'] = 0.0

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

def parseTranscriptInfo(transcript_dict):
    gene = transcript_dict['hgnc']
    biotype = transcript_dict['bioType']
    consequence = ','.join(transcript_dict['consequence'])
    transcriptID = transcript_dict['transcript']

    if 'proteinId' in transcript_dict:
        proteinID = transcript_dict['proteinId']
    else:
        proteinID = "-"

    if 'hgvsc' in transcript_dict:
        c_HGVS = transcript_dict['hgvsc']
    else:
        c_HGVS = "-"

    if 'hgvsp' in transcript_dict:
        p_HGVS = transcript_dict['hgvsp']
    else:
        p_HGVS = "-"

    if 'polyPhenScore' in transcript_dict:
        polyPhenScore = transcript_dict['polyPhenScore']
    else:
        polyPhenScore = '-'

    if 'polyPhenPred' in transcript_dict:
        polyPhenPred = transcript_dict['polyPhenPred']
    else:
        polyPhenPred = '-'

    if 'siftScore' in transcript_dict:
        siftScore = transcript_dict['siftScore']
    else:
        siftScore = '-'

    if 'siftPred' in transcript_dict:
        siftPred = transcript_dict['siftPrediction']
    else:
        siftPred = '-'

    info_string = f"{gene}|{transcriptID}|{c_HGVS}|{proteinID}|{p_HGVS}|{biotype}|{consequence}"

    return (info_string, polyPhenScore, polyPhenPred, siftScore, siftPred)

def parseNirvana(sample_id):
    main_dir = os.getcwd()
    output_dir = os.path.join(main_dir, "ParsedVariantReports")

    sample_output_dir = os.path.join(main_dir, f"{sample_id}_wgs")

    json_fname = get_json_variants_fname(sample_id, sample_output_dir)
    output_fname = get_output_csv_file(sample_id, output_dir)

    log_fname = os.path.join(sample_output_dir, f"{sample_id}.variant_parsing.log")
    sys.stdout.write(f"Parsing {json_fname}\n")

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
        with open(log_fname, 'w') as logfile:
            with gzip.open(json_fname, 'rt') as f:
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

                            if data['position_count'] % 100000 == 0:
                                sys.stdout.write(f"Processed {data['position_count']} positions\n")

                            if position_dict['filters'][0] == 'PASS':
                                out = defaultdict(None)
                                var_index = 0
                                allele_index = 1
                                samples_dict = position_dict['samples'][0]

                                out = parseBasicInfo(out, position_dict,samples_dict)

                                # Set a minimum coverage of 5X
                                try:
                                    depth = int(out.get('CovDepth'))
                                except:
                                    depth = -1
                                    logfile.write(f"No Coverage vale for variant at position: {out['Chr']}-{out['Pos']}-{out['Ref']}\n")

                                if depth >= 10:
                                    if 'variants' in position_dict:
                                        for var_dict in position_dict['variants']:
                                            out = parseBasicVarInfo(out, var_dict, samples_dict, var_index, allele_index)

                                            # Set VAF of 0.35 threshold
                                            if out.get('VAF') >= 0.35:
                                                # Set the 1% Allele Frequency cutoff here using the Controls
                                                if out.get('gnomAD_controlsAllAf') <= 0.005:
                                                    if 'transcripts' in var_dict:
                                                        info_strings = list()
                                                        polyPhenScore_list = list()
                                                        polyPhenPred_list = list()
                                                        siftScore_list = list()
                                                        siftPred_list = list()

                                                        for transcript_dict in var_dict['transcripts']:
                                                            (info_string, polyPhenScore, polyPhenPred, siftScore, siftPred) = parseTranscriptInfo(out, transcript_dict)

                                                            info_strings.append(info_string)
                                                            polyPhenScore_list.append(polyPhenScore)
                                                            polyPhenPred_list.append(polyPhenPred)
                                                            siftScore_list.append(siftScore)
                                                            siftPred_list.append(siftPred)

                                                        out['Gene and Transcript Info'] = ";".join(info_strings)
                                                        out['polyPhenScore'] = ";".join(polyPhenScore_list)
                                                        out['polyPhenPred'] = ";".join(polyPhenPred_list)
                                                        out['siftScore'] = ";".join(siftScore_list)
                                                        out['siftPred'] = ";".join(siftPred_list)

                                                        writer.writerow(out)
                                                    else:
                                                        writer.writerow(out)

                                            var_index += 1
                                            allele_index += 1
                                else:
                                    # sys.stderr.write(f"Position in Annotated JSON with no Variant entries: {out['Chr']}-{out['Pos']}-{out['Ref']}\n")
                                    logfile.write(f"Position in Annotated JSON with no Variant entries: {out['Chr']}-{out['Pos']}-{out['Ref']}\n")



                        if is_gene_line:
                            data['genes'].append(trimmed_line.rstrip(','))
                            data['gene_count'] += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples', help="Input file for samples")
    parser.add_argument('-t', '--threads', help="Number of Threads", default=1)

    args = parser.parse_args()
    args.logLevel = "INFO"

    pool = Pool(int(args.threads))
    commands = list()

    main_dir = os.getcwd()
    output_dir = os.path.join(main_dir, "ParsedVariantReports")

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    samples = list()
    with open(args.samples, 'r') as input:
        reader = csv.DictReader(input, delimiter=',', quotechar='"')
        for row in reader:
            samples.append(row['sampleID'])

    results = pool.map(parseNirvana, samples)
