'''
Name: Jairo Navarro
Date: June 27, 2015

This program takes in a dictionary and converts a genomic
sequence from one organism and codon optimizes it for another
organism

Usage: 

> python3 optimizer.py Dictionary1 Dictionary2 Sequence
> Enter genomic sequence file: sequence.fa

Output:
Codon optimized sequence for dictionary2 in a text file

Exit Codes:
  0: No errors
  1: Two dictionaries not given to program
  2: No file extension on the dictionary files
  3: File does not exist
  4: More than 3 parameters after the program name
'''
class optimizer:

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

    def __init__ (self, refDic, convDic):
        '''
        Containers for the reference dictionary
                           converstion dictionary
        '''
        # Save the file names
        self.refDic = refDic
        self.convDic = convDic
        #CodonBiasDicionaries
        self.CBDcon = {}
        self.CBDref = {}


    def OpenFile(self, seqFile):
        '''Open the sequence file and ensure that it exists.'''
        try:
           return open(seqFile,'r')
        except FileNotFoundError:
           print("\n%s does not exist.\n" % seqFile)
           sys.exit(3)

    def CodonToAmino(self, codon):
        # Change DNA to RNA and return and AA
        return self.dnaCodonTable[codon]

    '''    Codon | AA | Freq | Count   '''
    def parseTable(self, line):
        # Create a list by spliting by whitespace
        List = line.split()
        return List

    def MakeDict (self):
        # Open the conversion codon bias table
        ConFile = self.OpenFile(self.convDic)
        # Read line by line
        for line in ConFile:
            # List returned from parseTable function
            currList = self.parseTable(line)
            # Add codon info to the dictionary
            self.convertBD(currList)

        # Follow the same steps for the reference codon bias table
        RefFile = self.OpenFile(self.refDic)
        for line in RefFile:
            currList = self.parseTable(line)
            self.referenceBD(currList)
        # Sort the dictionaries 
        self.CBDref = self.sortDict(self.CBDref)
        self.CBDcon = self.sortDict(self.CBDcon)


    ''' convertBD and referenceBD are used to have the
        frequency and codon counts to be recognized as 
        floats and ints respectively.
    '''
    def convertBD (self,List):
        # Extract info from list and add to the dictionary
        codon = List[0]
        AA = List[1]
        freq = float(List[2])
        count = int(List[3])
        try:
            self.CBDcon[AA].append([codon,freq,count])
        except KeyError:
            self.CBDcon[AA] = [[codon,freq,count]]

    def referenceBD(self,List):
        # Extract info from list and add to the dictionary
        codon = List[0]
        AA = List[1]
        freq = float(List[2])
        count = int(List[3])
        try:
            self.CBDref[AA].append([codon,freq,count])
        except KeyError:
            self.CBDref[AA] = [[codon,freq,count]]

    def sortDict(self, dic):
        # Sort codons per AA by frequency
        # Highest frequency first

        from operator import itemgetter
        for entry,value in dic.items():
            sortedList = sorted(value, key=itemgetter(1), reverse=True)
            dic[entry] = sortedList
        return dic

    def recode(self,codon):
        # Get AA from the codon
        AA = self.CodonToAmino(codon)
        # Sorted codon list by frequency for the AA
        codonList = self.CBDref[AA]
        position = 0
        # Go through the list until codon is found
        for entry in codonList:
            if codon == entry[0]:
                # Break once position is found
                break
            position += 1
        # Return the codon that is in the same position for 
        # the other organism
        return self.CBDcon[AA][position][0]

def main():
    import fastaReader as fRead
    import sys
    #####################################################
    '''Replace with CommandLine class'''

    #First file given in command line after program name
    REF_ORG = sys.argv[1]
    #Second file given in the command line
    CONV_ORG = sys.argv[2]
    # FASTA file
    FASTA = sys.argv[3]
    OutFile = "Optimized-" + FASTA
    
    #####################################################
    FileOutput = open(OutFile,'w+')
    thisReader = fRead.FastAreader(FASTA)
    Recoder = optimizer(REF_ORG, CONV_ORG)

    # Make and sort codon bias tables by frequency
    Recoder.MakeDict()

    for head,seq in thisReader.readFasta():
        recodedSeq = ''

        for nuc in range(0, len(seq), 3):
            codon = seq[nuc:nuc+3]
            recodedSeq += Recoder.recode(codon)
        
        print('>' + head, file = FileOutput)
        FormattedSeq = ''
        count = 0
        for char in recodedSeq:
            if count == 70:
                FormattedSeq += "\n"
                count = 0
            FormattedSeq += char
            count += 1
        print(FormattedSeq, file = FileOutput)
        

if __name__ == "__main__":
    main()

