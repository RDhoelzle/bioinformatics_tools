#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = sys.argv[1] 
identifier = sys.argv[2]
result_file = identifier + ".fa"
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as output_file:
  for seq in fasta_sequences:
    if identifier in seq.id :
      SeqIO.write([seq], output_file, "fasta")
