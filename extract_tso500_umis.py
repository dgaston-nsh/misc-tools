#!/usr/bin/env python


# This script is used to extract the first 7 nucleotides from reads in de-multiplexed
# FastQ files and add them to the header in the 8th field. This is for cases where a
# sequencing centre ran a standard de-multiplexing workflow and now the proper BCL
# workflow for TSO500 data.


# Standard packages
import os
import glob
import HTSeq
import shutil

def processFastQ(filename, outfilename):
    with HTSeq.FastqReader(filename) as fq:
        with open(outfilename, 'w') as fq_out:
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
