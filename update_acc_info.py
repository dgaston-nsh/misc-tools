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

    case_data = defaultdict(lambda: dict)

    with open(args.update, 'r') as input:
        reader = csv.DictReader(input, delimiter=',')
        for row in reader:
            if args.type == 'WGS':
                
