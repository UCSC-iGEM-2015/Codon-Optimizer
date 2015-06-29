#!/user/bin/env python3
# Author: Jairo Navarro
# Date: June 26, 2015

'''
This program reads in a FASTA file(.fa) from the command line and 
creates a dictionary for all of the codons in an open reading
frame (ORF). The dictionary will then be written to a file with the 
same name as the input FASTA file with a .dic extension. 

Usage:
python3 read_genome.py Example.fa 

Output: 
Example.dic

'''

class read_genome :
    
    def __init__ (self, fileName = ''):
        
        #Saves the attribute file name
        self.fileName = fileName

    def OpenFile (self):

        # Opens the specified file
        if self.fileName is '':
            return sys.stdin
        else:
            return open(self.fileName)

    def ReadFile (self):
        header = ''
        sequence = ''
        
        with self.OpenFile() as fileH:
            # initialize return containers
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
 
            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence
