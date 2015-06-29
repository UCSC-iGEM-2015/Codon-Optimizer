#!/usr/bin/env python3
#Name: Jairo Navarro (jnavarr5)
#Group Members: 'None'

'''
This program has creates 3 different dictionaries. 
One gets the amount of codons that are present in the sequence provided (nucComposition).
The second gets the different codons that are present in the sequence (codonComposition).
The third one gets the total number of each amino acid that is encoded by the sequence (aaComposition)

The input that this program expects is a string that is longer than 3 characters in length 
The program does not print anything, but instead returns 3 different dictionaries that can be used in another program.
'''

class NucParams:

    rnaCodonTable = {

    # RNA codon table

    # U

    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', # UxU

    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C', # UxC

    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-', # UxA

    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W', # UxG

    # C

    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', # CxU

    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', # CxC

    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', # CxA

    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', # CxG

    # A

    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', # AxU

    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', # AxC

    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', # AxA

    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', # AxG

    # G

    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', # GxU

    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', # GxC
    
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', # GxA

    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G' # GxG

    }

    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    

    def __init__ (self):
        
        #Dictionaries for the genome
        self.nucDict = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0}
        self.codonDict = {}
        self.aaCompDict = {}

        #These next two for loops creates keys for the codons that were not in the input sequence.

        #Loops through the RNA Codon Table 
        for key in NucParams.rnaCodonTable.keys():
            #Gives the missing key the value 0.
            self.codonDict[key] = 0
            self.aaCompDict[NucParams.rnaCodonTable[key]] = 0


    def addSequence (self, thisSequence):
        import re
        #Creates a string all of the whitespace removed and all characters will be uppercase
        seq = thisSequence.upper().replace('\t','').replace(' ', '')
        #Creating a lis of all codons in the sequence(Splits at every 3 char)
        codonList = re.findall('...',seq)
        nucList = re.findall('.',seq)
        #Counter for removing invalid codons
        count = 0
        #Loops through the list to remove unwanted codons
        for codon in codonList:
            for base in codon:
                if base not in self.nucDict.keys():
                    count +=1
            if count >0:
                codonList.remove(codon)
                count = 0
                
        #This for loop with then go through the list and add values to the genome dictionary 
        for codon in codonList:
            if codon in NucParams.rnaCodonTable.keys():
                self.codonDict[codon] += 1
                self.aaCompDict[NucParams.rnaCodonTable[codon]] += 1
               
            elif codon in NucParams.dnaCodonTable.keys():
                self.codonDict[codon.replace('T','U')] += 1  
                self.aaCompDict[NucParams.rnaCodonTable[codon.replace('T','U')]] += 1
        #This for loop adds keys and values to the nucleotides dictionary
        for base in nucList:
            if base in self.nucDict.keys():
                self.nucDict[base] += 1

    def aaComposition(self):
    
        return self.aaCompDict

    def nucComposition(self):

        return self.nucDict  

    def codonComposition(self):

        return self.codonDict

    def nucCount(self,nuc):

        count = 0
        #Loops through the Nucleotide dictionary and gets the total amount of nucleotides
        for key in self.nucDict:
            count += self.nucDict[key]
        return self.nucDict[nuc]/count
    

 

 

class FastAreader :

    '''

    Class to provide reading of a file containing one or more FASTA

    formatted sequences:

    object instantiation:

    FastAreader(<file name>):

 

    object attributes:

    fname: the initial file name

 

    methods:

    readFasta() : returns header and sequence as strings.

    Author: David Bernick

    Date: April 19, 2013

    '''

    def __init__ (self, fname):

        '''contructor: saves attribute fname '''

        self.fname = fname

 

    def readFasta (self):

        '''

        using filename given in init, returns each included FastA record

        as 2 strings - header and sequence.

        whitespace is removed, no adjustment is made to sequence contents.

        The initial '>' is removed from the header.

        '''

        header = ''

        sequence = ''

        

        with open(self.fname) as fileH:

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

 

# presumed object instantiation and example usage

# myReader = FastAreader ('testTiny.fa');

# for head, seq in myReader.readFasta() :

#     print (head,seq)

