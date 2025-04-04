#!/usr/bin/env python


# This script is used to extract the first 7 nucleotides from reads in de-multiplexed
# FastQ files and add them to the header in the 8th field. This is for cases where a
# sequencing centre ran a standard de-multiplexing workflow and now the proper BCL
# workflow for TSO500 data.


# Standard packages
import os
import gzip
import glob
import HTSeq
import shutil

def extractUMIs_paired_end(fastq1, fastq2, outfilename1, outfilename2):
    with gzip.open(fastq1, 'rt') as fq1:
        with gzip.open(fastq2, 'rt') as fq2:
            with gzip.open(outfilename1, 'wt') as out1:
                with gzip.open(outfilename2, 'wt') as out2:
                    lines_fq1 = islice(fq1, 4)
                    lines_fq2 = islice(fq2, 4)




def processFastQ(fastq1, fastq2, outfilename1, outfilename2):
    with HTSeq.FastqReader(fastq1) as fq1:
        with open(outfilename1, 'w') as fq_out2:
            for read in fq:
                umi_record = read[:7]
                umi_seq = umi_record.read.decode()

                trimmed_record = read[7:]
                read_id_sections = read.name.split(' ')

                new_read_id = f"{read_id_sections[0]}:{umi_seq} {read_id_sections[1]}"
                trimmed_record.name = new_read_id

                trimmed_record.write_to_fastq_file(fq_out)


if __name__ == "__main__":
    cwd = os.getcwd()
