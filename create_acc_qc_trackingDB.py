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
    parser.add_argument('-d', '--db', help="File path/name for SQLite3 DB")

    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write(f"Creating SQLite DB at {args.db}\n")
    db_connection = sqlite3.connect(f"{args.db}")
    db_cursor = db_connection.cursor()

    db_cursor.execute("CREATE TABLE wgs_tumour(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    db_cursor.execute("CREATE TABLE wgs_normal(acc_id, external_id, cohort, sequencing_centre, coverage, status)")
    db_cursor.execute("CREATE TABLE transcriptome(acc_id, external_id, cohort, sequencing_centre, million_reads, status)")
