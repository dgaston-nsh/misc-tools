#!/usr/bin/env python

# Standard packages
import os
import csv
import sys
import sqlite3
import argparse

VERSION = "1.0"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--db', help="File path/name of case-level QC tracking database")
    parser.add_argument('-i', '--input', help="File with initial case-level QC tracking data")

    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write(f"Creating SQLite DB at {args.db}\n")
    db_connection = sqlite3.connect(f"{args.db}")
    db_cursor = db_connection.cursor()

    rna_inserts = list()
    norm_wgs_inserts = list()
    tumour_wgs_inserts = list()

    #("CREATE TABLE wgs_tumor(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    #("CREATE TABLE wgs_normal(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    #("CREATE TABLE transcriptome(acc_id, external_id, cohort, sequencing_centre, million_reads, status)")

    sys.stdout.write(f"Loading case-level QC tracking database with info from {args.input}\n")
    with open(args.input, 'r') as input:
        reader = csv.DictReader(input, delimiter=',')
        for row in reader:
            if row['# Million RNA Read Pairs'] != 'N/A' and row['# Million RNA Read Pairs'] != '':
                #print(row['# Million RNA Read Pairs'])
                rna_reads = float(row['# Million RNA Read Pairs'])
                if rna_reads > 100.0:
                    status = 'Tier A'
                else:
                    status = 'Tier B'
                rna_inserts.append((row['ACC-ID'], row['ExternalID'], row['Cohort'], row['Sequencing Centre'], rna_reads, status))

            if row['Normal WGS Coverage (30X)'] != 'N/A' and row['Normal WGS Coverage (30X)'] != '':
                coverage = float(row['Normal WGS Coverage (30X)'])
                if coverage >= 30.0:
                    status = 'Tier A'
                else:
                    status = 'Incomplete - Needs Top Up'
                norm_wgs_inserts.append((row['ACC-ID'], row['ExternalID'], row['Cohort'], row['Sequencing Centre'], coverage, status))

            if row['Tumour WGS Coverage (80X/30X)'] != 'N/A' and row['Tumour WGS Coverage (80X/30X)'] != '' and row['Tumour WGS Coverage (80X/30X)'] != 'pending':
                coverage = float(row['Tumour WGS Coverage (80X/30X)'])
                if coverage >= 80.0:
                    status = 'Tier A'
                elif coverage >= 30.0:
                    status = 'Tier B'
                else:
                    status = 'FAIL'
                tumour_wgs_inserts.append((row['ACC-ID'], row['ExternalID'], row['Cohort'], row['Sequencing Centre'], coverage, status))

    db_cursor.executemany("INSERT INTO transcriptome VALUES(?, ?, ?, ?, ?, ?)", rna_inserts)
    db_connection.commit()

    db_cursor.executemany("INSERT INTO wgs_normal VALUES(?, ?, ?, ?, ?, ?)", norm_wgs_inserts)
    db_connection.commit()

    db_cursor.executemany("INSERT INTO wgs_tumour VALUES(?, ?, ?, ?, ?, ?)", tumour_wgs_inserts)
    db_connection.commit()
