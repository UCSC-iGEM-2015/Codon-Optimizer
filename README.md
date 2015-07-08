# Codon-Optimizer
Codon optimizer written in Python for the University of California, Santa Cruz International genetically engineered machines competition (iGEM). 

To create a codon bias table, use percentageCount.py

  Usage: python3 percentageCount.py FASTA.fa
  
  All fasta files must be gene sequences with AUG and and stop codon for each gene sequence
  You need to run this twice to generate your source and target genome. 

To codon optimize a nucleic sequence, optimizer.py is used.

  Usage: python3 optimizer.py Reference.dic Conversion.dic
