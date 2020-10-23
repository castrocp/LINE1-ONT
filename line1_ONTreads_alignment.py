#!/usr/bin/env python3

###########################
#  Run this as:
# 'script_name.py L1.fasta long_reads_input_file.fasta outfile.txt'
#
##########################

import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

L1_fasta = sys.argv[1]
read_fasta = sys.argv[2]
outfile = sys.argv[3]

############################################
# Store L1 sequence and reverse complement #
############################################
L1 = SeqIO.read(L1_fasta, "fasta")
L1_seq = L1.seq
rev_comp_L1seq = L1_seq.reverse_complement()

##################
# Create Aligner #
##################

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 5
aligner.mismatch_score = -5
aligner.open_gap_score = -10
aligner.extend_gap_score = -5

##############################
# Align reads to L1 sequence #
##############################
with open(outfile, "w") as results_out:
    i = 1 # Keep track of how many reads have been processed
    
    for seq in SeqIO.parse(read_fasta, "fasta"):
        if i % 50 == 0:
            print(str(i) + " reads processed")
        results_out.write(seq.id + "\n")

        # Use Bio.Align package to do pairwise alignments

        forward_alignments = aligner.align(seq.seq[:80], L1_seq)
        if not any(a.score >= 100 for a in forward_alignments):
            results_out.write("No alignments with score >= 100 on the forward strand" + "\n")
        else:
            if ((forward_alignments[0].aligned)[0][0][0]) == 0:
                results_out.write("Forward alignment score: " + "\t" + str(forward_alignments[0].score) + "\n")
                results_out.write(str(forward_alignments[0].aligned) + "\n")
                results_out.write("Forward alignment starting L1HS base:" + "\t" + str((forward_alignments[0].aligned)[1][0][0] + 1) + "\n")
            else:
                results_out.write("Forward strand alignment does not begin at start of read" + "\n")

        rev_comp_alignments = aligner.align(seq.seq[:80], rev_comp_L1seq)
        if not any(a.score >= 100 for a in rev_comp_alignments):
            results_out.write("No alignments with score >= 100 on the reverse strand" + "\n")
            results_out.write("\n")
        else:
            if ((rev_comp_alignments[0].aligned)[0][0][0]) == 0:
                results_out.write("Reverse strand alignment score: " + "\t" + str(rev_comp_alignments[0].score) + "\n")
                results_out.write(str(rev_comp_alignments[0].aligned) + "\n")
                results_out.write("Reverse alignment starting L1HS base:" + "\t" + str( 6059 - ((rev_comp_alignments[0].aligned)[1][0][0])) + "\n")
                results_out.write("\n")
            else:
                results_out.write("Reverse strand alignment does not begin at start of read" + "\n")
                results_out.write("\n")
                
        i += 1
  
