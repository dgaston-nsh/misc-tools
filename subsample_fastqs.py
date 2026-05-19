#!/usr/bin/env python

import os
import gzip
import argparse
from multiprocessing import Pool

def subset_fastq(args):
    """
    Subsets a single FASTQ file to a specific number of reads.
    'args' is a tuple of (input_path, output_path, num_reads)
    """
    input_file, output_file, target_reads = args

    # Target reads in millions; convert to actual count
    max_records = int(target_reads * 1_000_000)

    print(f"Processing: {os.path.basename(input_file)}...")

    try:
        # Open both files in streaming mode (binary)
        with gzip.open(input_file, 'rt') as f_in, \
             gzip.open(output_file, 'wt', compresslevel=6) as f_out:

            count = 0
            # Read 4 lines at a time (one FASTQ record)
            while count < max_records:
                line1 = f_in.readline()
                if not line1: break # End of file
                line2 = f_in.readline()
                line3 = f_in.readline()
                line4 = f_in.readline()

                # Write the 4-line record to the output
                f_out.write(line1 + line2 + line3 + line4)

                count += 1

        print(f"Finished: {os.path.basename(input_file)} ({count} records)")
    except Exception as e:
        print(f"Error processing {input_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Parallel FASTQ subsetting tool.")
    parser.add_argument("-i", "--input", required=True, help="Input directory containing .fastq.gz files")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-n", "--num_reads", type=float, default=1.0, help="Number of million reads to keep (e.g., 0.5 or 2)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of parallel threads")

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Gather all .fastq.gz files
    tasks = []
    for filename in os.listdir(args.input):
        if filename.endswith(".fastq.gz"):
            in_path = os.path.join(args.input, filename)
            out_path = os.path.join(args.output, filename)
            tasks.append((in_path, out_path, args.num_reads))

    # Initialize a Pool of workers for parallel execution
    print(f"Starting subsetting with {args.threads} threads...")
    with Pool(processes=args.threads) as pool:
        pool.map(subset_fastq, tasks)

    print("All tasks completed.")

if __name__ == "__main__":
    main()
