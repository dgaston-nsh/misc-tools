#!/usr/bin/env python

import csv
import sys
import sqlite3
import argparse

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help="Output file name")
    parser.add_argument('-d', '--db', help="File path/name of case-level QC tracking database")

    args = parser.parse_args()
    args.logLevel = "INFO"

    db_connection = sqlite3.connect(f"{args.db}")
    db_cursor = db_connection.cursor()

    #("CREATE TABLE wgs_tumor(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    #("CREATE TABLE wgs_normal(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    #("CREATE TABLE transcriptome(acc_id, external_id, cohort, sequencing_centre, million_reads, status)")

    db_cursor.execute("SELECT * FROM wgs_tumour")
    tumour_results = db_cursor.fetchall()

    db_cursor.execute("SELECT * FROM wgs_normal")
    normal_results = db_cursor.fetchall()

    db_cursor.execute("SELECT * FROM transcriptome")
    transcriptome_results = db_cursor.fetchall()

    case_data = defaultdict(lambda: defaultdict())

    for res in tumour_results:
        case_data[res[0]]['ExternalID'] = res[1]
        case_data[res[0]]['Cohort'] = res[2]
        case_data[res[0]]['Sequencing Centre'] = res[3]
        case_data[res[0]]['Tumour Coverage'] = res[4]
        case_data[res[0]]['Tumour Status'] = res[5]

    for res in normal_results:
        case_data[res[0]]['Normal Coverage'] = res[4]
        case_data[res[0]]['Normal Status'] = res[5]

        if case_data[res[0]].get('ExternalID') and res[1]:
            if case_data[res[0]]['ExternalID'] != res[1]:
                sys.stderr.write(f"{res[0]} Normal has mismatches in External IDs ({case_data[res[0]]['ExternalID']} and {res[1]})\n")
                case_data[res[0]]['ExternalID'] = res[1]
        else:
            if res[1] != '':
                case_data[res[0]]['ExternalID'] = res[1]

        if case_data[res[0]].get('Cohort') and res[2]:
            if case_data[res[0]]['Cohort'] != res[2]:
                sys.stderr.write(f"{res[0]} Normal has mismatches in Cohort designation ({case_data[res[0]]['Cohort']} and {res[2]})\n")
                case_data[res[0]]['Cohort'] = res[2]
        else:
            if res[2] != '':
                case_data[res[0]]['Cohort'] = res[2]

        if case_data[res[0]].get('Sequencing Centre') and res[3]:
            if case_data[res[0]]['Sequencing Centre'] != res[3]:
                sys.stderr.write(f"{res[0]} Normal has mismatches in Sequencing Centre ({case_data[res[0]]['Sequencing Centre']} and {res[3]})\n")
                case_data[res[0]]['Sequencing Centre'] = res[3]
        else:
            if res[3] != '':
                case_data[res[0]]['Sequencing Centre'] = res[3]

    for res in transcriptome_results:
        case_data[res[0]]['Millions of Reads'] = res[4]
        case_data[res[0]]['RNA Status'] = res[5]

        if case_data[res[0]].get('ExternalID') and res[1]:
            if case_data[res[0]]['ExternalID'] != res[1]:
                sys.stderr.write(f"{res[0]} Transcriptome has mismatches in External IDs ({case_data[res[0]]['ExternalID']} and {res[1]})\n")
                case_data[res[0]]['ExternalID'] = res[1]
        else:
            if res[1] != '':
                case_data[res[0]]['ExternalID'] = res[1]

        if case_data[res[0]].get('Cohort') and res[2]:
            if case_data[res[0]]['Cohort'] != res[2]:
                sys.stderr.write(f"{res[0]} Transcriptome has mismatches in Cohort designation ({case_data[res[0]]['Cohort']} and {res[2]})\n")
                case_data[res[0]]['Cohort'] = res[2]
        else:
            if res[2] != '':
                case_data[res[0]]['Cohort'] = res[2]

        if case_data[res[0]].get('Sequencing Centre') and res[3]:
            if case_data[res[0]]['Sequencing Centre'] != res[3]:
                sys.stderr.write(f"{res[0]} Transcriptome has mismatches in Sequencing Centre ({case_data[res[0]]['Sequencing Centre']} and {res[3]})\n")
                case_data[res[0]]['Sequencing Centre'] = res[3]
        else:
            if res[3] != '':
                case_data[res[0]]['Sequencing Centre'] = res[3]

    with open(args.output, 'w') as output:
        output.write("ACC ID, External ID, Cohort, Normal Coverage (30X), Normal Status, Tumour Coverage (80X/30X), Tumour Status, Millions of Read Pairs, RNA Status, Sequencing Centre\n")
        for case in case_data:
            output.write(f"{case}, {case_data[case].get('ExternalID')}, {case_data[case].get('Cohort')}, {case_data[case].get('Normal Coverage')}, {case_data[case].get('Normal Status')}, {case_data[case].get('Tumour Coverage')}, {case_data[case].get('Tumour Status')}, {case_data[case].get('Millions of Reads')}, {case_data[case].get('RNA Status')}, {case_data[case].get('Sequencing Centre')}\n")
