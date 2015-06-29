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
        self.refDic = refDic
        self.convDic = convDic
        self.dicts_from_ref = {}
        self.dicts_from_conv = {}
        self.codonList = []

    def ReadDic (self):
        '''
        Read the dictionary files and save to the containers
        '''
        try:
            for line in open(refDic,'r'):
                self.dicts_from_ref = eval(line)
        except FileNotFoundError:
            print("\n%s does not exist.\n" % refDic)
            sys.exit(3)
        try:
            for line in open(convDic,'r'):
                self.dicts_from_conv = eval(line)
        except FileNotFoundError:
            print("\n%s does not exist.\n" % convDic)
            sys.exit(3)
        #Returns a tuple of dictionaries
        return self.dicts_from_ref, self.dicts_from_conv

    def OpenSeqFile(self, seqFile):
        '''Open the sequence file and ensure that it exists.'''
        try:
           return open(seqFile,'r')
        except FileNotFoundError:
           print("\n%s does not exist.\n" % seqFile)
           sys.exit(3)
    def CodonToAmino(self, codon):
        return self.rnaCodonTable[codon.replace('T','U')]

    def GetValue(self, dic, AA,codon):
        from operator import itemgetter
        print(type(AA))
        print(type(dic))
        print(dic[AA])
        List = dic[AA]
        for count,entry in enumerate(List):
            if AA in entry:
                print(count)
                
        Value = dic[AA]
        #print(Value)
        #Value = Value[1]
        #Value = Value.split(']')
        #Value = sorted(eval("["+Value[0]+"]"), key=itemgetter(1),reverse=True)
        print(type(Value))
        for entry in Value:
            print(entry)
            print(type(entry))
            if codon in entry:
                print(codon)
                print(entry[codon])
                self.codonList.append((AA,entry[codon]))
        return self.codonList

    def AnalyzeRef (self,seq):
        '''
        Analyze the reference dictionary and return the AA and 
        optimal codon frequency within the reference genome.

        return value is a tuple.
        '''
        protein = ''
        seq = seq.replace('T','U')
        for NucAcid in range (0,len(seq),3):
            codon = seq[NucAcid : NucAcid +3]
            print(codon)
            print(self.refDic)
            protein = self.CodonToAmino(codon)
            print(protein)
            self.GetValue(self.dicts_from_ref,protein,codon)
        return self.codonList  

    def Translate (self,AAfreq):
        '''
        Takes in a list of tuple with the first value being the AA and the
        second value being its frequency in the reference genome
        '''
        finalSeq = ''
        for entry in AAfreq:
            AA = entry[0]
            position = entry[1]
            finalSeq += self.dicts_from_conv[AA][position]
        return finalSeq

    def AATranslate(self, AA, position):
        pass

import sys
seqFile = ''
refDic = ''
convDic = ''

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
refDic = sys.argv[1]
#Second file given in the command line
convDic = sys.argv[2]

#Make sure the dictionaries have a .dic extension
try:
    refDic.split('.')[1]
    convDic.split('.')[1]
except IndexError:
    print("\nError: Dictionaries must end in an extension.\n")
    print(usage)
    sys.exit(2)

#Create an object for the two dictionaries
mydic = read_dic(refDic,convDic)

#ref and conv are both dictionaries of dictionaries
#First key is the AA second key is codon
ref = mydic.ReadDic()[0]
conv = mydic.ReadDic()[1]
#Access a dictionary file name
print(mydic.convDic)
#Access a AA within a dictionary
print(ref["A"])

protein = mydic.AnalyzeRef("ATGATCTATAAGTAA")
print(protein)
''' 
Need to make sure that the sequence file exists
Add computations for the dictionary.
'''
