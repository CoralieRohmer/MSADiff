#!/usr/bin/env python3
import sys,re,os
from Bio import SeqIO

try:
  file_clustal=sys.argv[sys.argv.index("-in") + 1]
except:
  pass
  print("No file clustal format. -in file")
  sys.exit()

try:
  output_fasta=sys.argv[sys.argv.index("-out") + 1]
except:
  pass
  print("No output. -out file")
  sys.exit()

records = SeqIO.parse(file_clustal, "clustal")
count = SeqIO.write(records, output_fasta, "fasta")
