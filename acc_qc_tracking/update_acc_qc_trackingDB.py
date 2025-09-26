#!/usr/bin/env python

import csv
import sys
import sqlite3
import argparse

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--update', help="Update file to add")
    parser.add_argument('-t', '--type', help="Type of data update (WGS or WTS)")
    parser.add_argument('-d', '--db', help="File path/name of case-level QC tracking database")

    args = parser.parse_args()
    args.logLevel = "INFO"

    db_connection = sqlite3.connect(f"{args.db}")
    db_cursor = db_connection.cursor()

    #("CREATE TABLE wgs_tumor(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    #("CREATE TABLE wgs_normal(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    #("CREATE TABLE transcriptome(acc_id, external_id, cohort, sequencing_centre, million_reads, status)")

    rna_inserts = list()
    norm_wgs_inserts = list()
    tumour_wgs_inserts = list()

    with open(args.update, 'r') as input:
        reader = csv.DictReader(input, delimiter=',')
        for row in reader:
            query_params = (row['ACC-ID'],)
            if args.type == 'WGS':
                if row['Sample type'] == 'Tumor':
                    db_cursor.execute("SELECT cohort FROM wgs_tumour WHERE acc_id = ?", query_params)
                    results = db_cursor.fetchall()
                    if len(results) > 0:
                        sys.stderr.write(f"Data for {row['ACC-ID']} already exists in database, perform separate update\n")
                    else:
                        coverage = float(row['BWA_Aligned_dedup_cov'])
                        if coverage >= 80.0:
                            status = 'Tier A'
                        elif coverage >= 30.0:
                            status = 'Tier B'
                        else:
                            status = 'FAIL'

                        tumour_wgs_inserts.append((row['ACC-ID'], '', row['Cohort'], row['Sequencing Centre'], coverage, status))
                elif row['Sample type'] == 'Normal':
                    db_cursor.execute("SELECT cohort FROM wgs_normal WHERE acc_id = ?", query_params)
                    results = db_cursor.fetchall()
                    if len(results) > 0:
                        sys.stderr.write(f"Data for {row['ACC-ID']} already exists in database, perform separate update\n")
                    else:
                        coverage = float(row['BWA_Aligned_dedup_cov'])
                        if coverage >= 30.0:
                            status = 'Tier A'
                        else:
                            status = 'FAIL'

                        norm_wgs_inserts.append((row['ACC-ID'], '', row['Cohort'], row['Sequencing Centre'], coverage, status))
                else:
                    sys.stderr.write(f"Unknown sample type for {row['ACC-ID']}. Skipping...\n")

            elif args.type == 'Transcriptome':
                if row['Type (Tumor/Normal)'] == 'Tumor':
                    db_cursor.execute("SELECT cohort FROM transcriptome WHERE acc_id = ?", query_params)
                    results = db_cursor.fetchall()
                    if len(results) > 0:
                        sys.stderr.write(f"Data for {row['ACC-ID']} already exists in database, perform separate update\n")
                    else:
                        num_reads = float(row['# Million Reads'])
                        if num_reads >= 100.0:
                            status = 'Tier A'
                        else:
                            status = 'Tier B'

                        rna_inserts.append((row['ACC-ID'], '', row['Cohort'], row['Sequencing Centre'], num_reads, status))
                else:
                    sys.stderr.write(f"Non-Tumour sample type {row['Type (Tumor/Normal)']} for sample {row['ACC-ID']}. Skipping...\n")

            else:
                sys.stderr.write(f"Unknown data type {args.type}. Exiting...\n")
                sys.exit()

    if args.type == 'WGS':
        db_cursor.executemany("INSERT INTO wgs_normal VALUES(?, ?, ?, ?, ?, ?)", norm_wgs_inserts)
        db_connection.commit()

        db_cursor.executemany("INSERT INTO wgs_tumour VALUES(?, ?, ?, ?, ?, ?)", tumour_wgs_inserts)
        db_connection.commit()
    elif args.type == 'Transcriptome':
        db_cursor.executemany("INSERT INTO transcriptome VALUES(?, ?, ?, ?, ?, ?)", rna_inserts)
        db_connection.commit()
    else:
        sys.stderr.write(f"Unknown data type {args.type}. Exiting...\n")
        sys.exit()
