#!/usr/bin/env python3

import sys
import os
import subprocess
import pysam

def main():


project_dir = "/data/data_repo/castrocp/LINE1" 
ref = project_dir + "/data/hg38.fa"
ref_idx = project_dir + "/data/hg38.mmi"
reads_fastq = project_dir + "/data/aba133_e464c2c4703faa4988dfffa0e8beb1294afa9c0b_0.fastq" # make this a user argument
sam = project_dir + "/results/" + os.path.splitext(os.path.basename(reads_fastq))[0] + ".sam"
sorted_sam = project_dir + "/results/" + os.path.splitext(os.path.basename(sam))[0] + ".sorted.sam"
bam = project_dir + "/results/" + os.path.splitext(os.path.basename(sam))[0] + ".bam"


def create_index():
    index_cmd = ["minimap2", "-x", "map-ont", "-d", ref_idx, ref]
    subprocess.run(index_cmd, check=True)

def align():
    minimap_cmd = ["minimap2", "-ax", "map-ont", ref_idx, reads_fastq, "-o", out_sam]
    subprocess.run(minimap_cmd, check=True)

def compress():
    # Sort aligned reads
    pysam.sort("-o", sorted_sam, sam, catch_stdout=False)
    # Compress SAM into BAM format
    pysam.view("-b", "-o", bam, sorted_sam, catch_stdout=False)


def main()
    do the thing



if __name__ == "__main__":
    main()
