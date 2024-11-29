#!/usr/bin/env python

import os
import csv
import sys
import gzip
import json
import shutil
import fnmatch
import argparse

import nirvana

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help="Output directory")
    args = parser.parse_args()
    args.logLevel = "INFO"

    header = "SampleID\tVarID\tGene\tg.HGVS\tVarType\tTranscript\tBiotype\tCodons \
    \tAminoAcids\tcDNAPos\taaPos\tcdsPos\tExons\tIntrons\tc.HGVS\tp.HGVS\tdbSNP \
    \tFILTERS\tVAFS\talleleDepths \
    \tConsequence\tCancerResidue\tCancerNumSamples\tCancerNumAltAASamples\tCanceqValue \
    \tClinvarIDs\tClinvarSigs\tCOSMIC_IDs \
    \tgnomadAF\tTopmedAF\tPolyPhen\tSIFT\tDANN\tGERP\tREVEL\n"

    for root, dirs, files in os.walk("."):
        for nirvana_file in fnmatch.filter(files, "*.hard-filtered.vcf.annotated.json.gz"):

            elements = nirvana_file.split(".hard-filtered.")
            sample_name = elements[0]

            sys.stdout.write(f"Parsing JSON structure of file {nirvana_file}\n")
            data = nirvana.parseNirvana(os.path.join(root, nirvana_file).rstrip())
            sys.stdout.write("Processing and parsing position data\n")

            positions = 0
            variants = 0
            variants_written = 0
            not_pass = 0

            output_name = os.path.join(args.output, f"{sample_name}.filtered.pass.rare.variants.tsv")

            with open(output_name, 'w') as outfile:
                outfile.write(f"{header}")
                for position_dict in data['positions']:
                    positions += 1
                    if 'variants' in position_dict:
                        var_in_position = 0
                        if position_dict['filters'][0] == "PASS":
                            position_dict['samples'][var_in_position]
                            for var_dict in position_dict['variants']:
                                position_samples_dict = position_dict['samples'][var_in_position]
                                var_in_position += 1
                                if 'variantFrequencies' in position_samples_dict:
                                    variants += 1

                                    if 'transcripts' in var_dict:
                                        output_transcript_dicts = []
                                        for transcript_dict in var_dict['transcripts']:
                                            # Check if Transcript is in MANE + Clinical and Output
                                            if transcript_dict['transcript'].startswith("ENST"):
                                                germline_rare = True
                                                if 'gnomad' in var_dict:
                                                    if var_dict['gnomad']['allAf'] >= 0.005:
                                                        germline_rare = False
                                                if germline_rare:
                                                    variants_written += 1
                                                    varid = f"{position_dict['chromosome']}-{position_dict['position']}-{position_dict['refAllele']}-{var_dict['altAllele']}"

                                                    outfile.write(f"{data['samples'][0]}\t")

                                                    # Basic Variant Info
                                                    outfile.write(f"{varid}\t")
                                                    outfile.write(f"{transcript_dict['hgnc']}\t")
                                                    outfile.write(f"{var_dict['hgvsg']}\t")
                                                    outfile.write(f"{var_dict['variantType']}\t")
                                                    outfile.write(f"{transcript_dict['transcript']}\t")
                                                    outfile.write(f"{transcript_dict['bioType']}\t")
                                                    if 'codons' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['codons']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'aminoAcids' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['aminoAcids']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'cdnaPos' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['cdnaPos']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'proteinPos' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['proteinPos']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'cdsPos' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['cdsPos']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'exons' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['exons']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'introns' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['introns']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'hgvsc' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['hgvsc']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                    if 'hgvsp' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['hgvsp']}\t")
                                                    else:
                                                        outfile.write("-\t")

                                                    # outfile.write(f"{transcript_dict['isCanonical']}\t")

                                                    if 'dbsnp' in var_dict:
                                                        outfile.write(f"{','.join(var_dict['dbsnp'])}\t")
                                                    else:
                                                        outfile.write("-\t")

                                                    outfile.write(f"{','.join(position_dict['filters'])}\t")
                                                    outfile.write(f"{','.join([str(i) for i in position_samples_dict['variantFrequencies']])}\t")
                                                    outfile.write(f"{','.join([str(i) for i in position_samples_dict['alleleDepths']])}\t")

                                                    # Impact and Consequences
                                                    consequence = ",".join(transcript_dict['consequence'])
                                                    outfile.write(f"{consequence}\t")

                                                    if 'cancerHotspots' in transcript_dict:
                                                        outfile.write(f"{transcript_dict['cancerHotspots']['residue']}\t")
                                                        outfile.write(f"{transcript_dict['cancerHotspots']['numSamples']}\t")
                                                        outfile.write(f"{transcript_dict['cancerHotspots']['numAltAminoAcidSamples']}\t")
                                                        outfile.write(f"{transcript_dict['cancerHotspots']['qValue']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                        outfile.write("-\t")
                                                        outfile.write("-\t")
                                                        outfile.write("-\t")

                                                    if 'clinvar' in var_dict:
                                                        clinvar_ids = []
                                                        clinvar_sigs = []
                                                        for clinvar_dict in var_dict['clinvar']:
                                                            clinvar_ids.append(clinvar_dict['id'])
                                                            clinvar_sigs.append(','.join(clinvar_dict['significance']))
                                                        outfile.write(f"{','.join(clinvar_ids)}\t")
                                                        outfile.write(f"{';'.join(clinvar_sigs)}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                        outfile.write("-\t")

                                                    if 'cosmic' in var_dict:
                                                        cosmic_ids = []
                                                        sample_counts = []
                                                        for cosmic_dict in var_dict['cosmic']:
                                                            cosmic_ids.append(cosmic_dict['id'])
                                                            # sample_counts.append(cosmic_dict['sampleCount'])
                                                        outfile.write(f"{','.join(cosmic_ids)}\t")
                                                        # outfile.write(f"{','.join(sample_counts)}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                        # outfile.write("-\t")

                                                    if 'gnomad' in var_dict:
                                                        outfile.write(f"{var_dict['gnomad']['allAf']}\t")
                                                        # outfile.write(f"{var_dict['gnomad']['lowComplexityRegion']}\t")
                                                        # outfile.write(f"{var_dict['gnomad']['failedFilter']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                        # outfile.write("-\t")
                                                        # outfile.write("-\t")

                                                    if 'topmed' in var_dict:
                                                        outfile.write(f"{var_dict['topmed']['allAf']}\t")
                                                        # outfile.write(f"{var_dict['gnomad']['failedFilter']}\t")
                                                    else:
                                                        outfile.write("-\t")
                                                        # outfile.write("-\t")

                                                    # if 'gnomadExome' in var_dict:
                                                    #     outfile.write(f"{vardict_dict['gnomadExome']['']}\t")
                                                    # else:
                                                    #     outfile.write("-\t")

                                                    #Predictions
                                                    if 'polyPhenPrediction' in var_dict:
                                                        outfile.write(f"{transcript_dict['polyPhenPrediction']}\t")
                                                    else:
                                                        outfile.write("-\t")

                                                    if 'siftPrediction' in var_dict:
                                                        outfile.write(f"{transcript_dict['siftPrediction']}\t")
                                                    else:
                                                        outfile.write("-\t")

                                                    if 'dannScore' in var_dict:
                                                        outfile.write(f"{var_dict['dannScore']}\t")
                                                    else:
                                                        outfile.write("-\t")

                                                    if 'gerpScore' in var_dict:
                                                        outfile.write(f"{var_dict['gerpScore']}\t")
                                                    else:
                                                        outfile.write("-\t")

                                                    if 'revel' in var_dict:
                                                        outfile.write(f"{var_dict['revel']['score']}\n")
                                                    else:
                                                        outfile.write("-\n")

                                                    # outfile.write(f"{transcript_dict['']}\t")
                            var_in_position += 1
            sys.stdout.write(f"Processed {variants} variants across {positions} genomic positions and wrote {variants_written} PASS variants with transcript data and < 0.5% frequency in gnomad to file\n")
