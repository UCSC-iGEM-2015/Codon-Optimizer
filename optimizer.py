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
class read_dic:
    import sys
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
        # list of tuples
        self.codonList = []
        self.header = ''


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

    def AnalyzeRef (self,seqFile):
        '''
        Analyze the reference dictionary and return the AA and 
        optimal codon frequency within the reference genome.

        return value is a List of tuples.
        '''
        import re
        seq = ''
        # Open the Nucleotide sequence file
        seqFile = self.OpenFile(seqFile)
        # Read line by line
        for line in seqFile:
            # Save the header
            if line.startswith('>'):
                self.header = line.replace("\n","")
                continue
            # Remove all whitespace
            seq += ''.join(line.split())

        '''Seq now is the entire nucleotide sequence from the file'''
        # Read through the sequence, 3 characters at a time
        for NucAcid in range (0,len(seq),3):
            codon = seq[NucAcid : NucAcid +3]
            # Append to the codon list the tuple from GetValue
            # (codon,#)
            self.codonList.append(self.GetValue(self.CBDref,codon))
        return self.codonList   

    def GetValue(self,dic,codon):
        # Get AA from the codon
        AA = self.CodonToAmino(codon)
        # Sorted codon list by frequency for the AA
        codonList = dic[AA]
        freq = 0.0
        position = 0
        for entry in codonList:
            if codon == entry[0]:
                freq = entry[1]
                # return codon and the position in the codon list
                return (codon,position)
            position += 1

    def Translate (self):
        '''
        Takes in a list of tuple with the first value being the AA and the
        second value being its frequency in the reference genome
        '''
        # Final optimized sequence container
        finalSeq = ''
        # sort the converstion dictionary
        self.sortDict(self.CBDcon)
        count = 0
        # Read through the list of tuples
        # (codon,position)
        for entry in self.codonList:
            # Add a newline character to the sequence after 23 AA
            if count == 23:
                finalSeq += '\n'
                count = 0
            codon = entry[0]
            position = entry[1]
            # Add optimized sequence to the string
            finalSeq += self.optimize(codon,position)
            count += 1
        return finalSeq

    def optimize(self, codon, position):
        AA = self.CodonToAmino(codon)
        # Get the list of codons for the AA
        codonList = self.CBDcon[AA]
        # Return the codon with the closest frequency
        return codonList[position][0]

'''
  Main Program
'''

import sys
seqFile = ''

usage = "\nUsage: \npython3 %s Reference_dictionary Conversion_dictionary\n" % sys.argv[0]

usage2 = usage + "Optional: Sequence file after second dictionary.\n"

#Ensure that two dictionaries were given
if len(sys.argv) < 3:
    print("\nPlease enter two dictionary files.")
    print(usage)
    sys.exit(1)

if len(sys.argv) > 4:
    print(usage2)
    sys.exit(4)

if len(sys.argv) == 4:
    seqFile = sys.argv[3]
else: 
    seqFile = input("Enter the sequence file to be optimized:\n")


#First file given in command line after program name
refFile = sys.argv[1]
#Second file given in the command line
convFile = sys.argv[2]

# Give file names to the read_dic class
mydic = read_dic(refFile,convFile)

# Create dictionaries from the files
mydic.MakeDict()

# sort the Codon bias tables 
mydic.sortDict(mydic.CBDcon)
mydic.sortDict(mydic.CBDref)

# Analyze the sequence to be optimized
mydic.AnalyzeRef(seqFile)

# Optimized sequence saved as Translation
Translation = mydic.Translate()

# Write to file.
OUTFILE = "Optimized-"+seqFile
FileOutput = open(OUTFILE, 'w+')
print("New nucleic acid sequence saved as: \n%s" % OUTFILE)
print(mydic.header,file = FileOutput)
print(Translation,file = FileOutput)
