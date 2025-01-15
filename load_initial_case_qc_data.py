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

    sys.stdout.write(f"Loading case-level QC tracking database with info from {args.input}\n")
    with open(args.input, 'r') as input:
        reader = csv.DictReader(input, delimiter=',')
        for row in reader:
