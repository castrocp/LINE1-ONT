#!/usr/bin/env python3

import sys
import os
import subprocess
import pysam


project_dir = "/data/data_repo/castrocp/LINE1" 
ref = project_dir + "/data/hg38.fa"
ref_idx = project_dir + "/data/hg38.mmi"
reads_fastq = project_dir + "/data/aba133_e464c2c4703faa4988dfffa0e8beb1294afa9c0b_0.fastq"
sam = project_dir + "/results/" + os.path.splitext(os.path.basename(reads_fastq))[0] + ".sam"
sorted_sam = project_dir + "/results/" + os.path.splitext(os.path.basename(sam))[0] + ".sorted.sam"
bam = project_dir + "/results/" + os.path.splitext(os.path.basename(sam))[0] + ".bam"

print(ref +"\n" +  ref_idx + "\n" +  reads_fastq + "\n" + sam, sorted_sam, bam)
