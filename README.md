# Codon-Optimizer
Codon optimizer written in Python for the University of California, Santa Cruz International genetically engineered machines competition (iGEM). 

To create a codon bias table, percentageCount.py is used.
                       Usage:
       python3 percentageCount.py FASTA.fa
Where FASTA.fa is a FASTA file that contains the genome for the reference organism.


To codon optimize a nucleic sequence, optimizer.py is used.
                       Usage:
    python3 optimizer.py Reference.dic Conversion.dic
Where Reference.dic and Conversion.dic were created from the percentCount.py program. 
