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
        #CodonBiasDicionaries
        self.CBDcon = {}
        self.CBDref = {}
        self.codonList = []

    '''
    def ReadDic (self):
        #Read the dictionary files and save to the containers
        try:
            for line in open(self.refDic,'r'):
                self.CBDref = eval(line)
        except FileNotFoundError:
            print("\n%s does not exist.\n" % refDic)
            sys.exit(3)
        try:
            for line in open(self.convDic,'r'):
                self.CBDcon = eval(line)
        except FileNotFoundError:
            print("\n%s does not exist.\n" % convDic)
            sys.exit(3)
        #Returns a tuple of dictionaries
        return self.CBDref, self.CBDcon
    '''

    def OpenFile(self, seqFile):
        '''Open the sequence file and ensure that it exists.'''
        try:
           return open(seqFile,'r')
        except FileNotFoundError:
           print("\n%s does not exist.\n" % seqFile)
           sys.exit(3)

    def CodonToAmino(self, codon):
        # Change DNA to RNA and return and AA
        return self.rnaCodonTable[codon.replace('T','U')]

    '''    Codon | AA | Freq | Count   '''
    def parseTable(self, line):
        # Create a list by spliting by whitespace
        List = line.split()
        return List

    def MakeDict (self):
        try:
            
            ConFile = self.OpenFile(self.convDic)
            for line in ConFile:
                currList = self.parseTable(line)
                self.convertBD(currList)
            RefFile = self.OpenFile(self.refDic)
            for line in RefFile:
                currList = self.parseTable(line)
                self.referenceBD(currList)
        except KeyError:
            pass


    def convertBD (self,List):
        codon = List[0].replace("T","U")
        AA = List[1]
        freq = float(List[2])
        count = int(List[3])
        try:
            self.CBDcon[AA].append([codon,freq,count])
        except KeyError:
            self.CBDcon[AA] = [[codon,freq,count]]
    def referenceBD(self,List):
        codon = List[0].replace("T","U")
        AA = List[1]
        freq = float(List[2])
        count = int(List[3])
        try:
            self.CBDref[AA].append([codon,freq,count])
        except:
            self.CBDref[AA] = [[codon,freq,count]]

    def sortDict(self, dic):
        from operator import itemgetter
        for entry,value in dic.items():
            sortedList = sorted(value, key=itemgetter(1), reverse=True)
            dic[entry] = sortedList

    def AnalyzeRef (self,seqFile):
        '''
        Analyze the reference dictionary and return the AA and 
        optimal codon frequency within the reference genome.

        return value is a tuple.
        '''
        import re
        seq = ''
        seqFile = self.OpenFile(seqFile)
        for line in seqFile:
            if line.startswith('>'):
                continue
            seq += re.sub(r"[^\w\s]", "", line)
            seq = seq.replace("\n", '')
        seq = seq.replace('T','U')
        for NucAcid in range (0,len(seq),3):
            codon = seq[NucAcid : NucAcid +3]
            self.codonList.append(self.GetValue(self.CBDref,codon))
        return self.codonList   

    def GetValue(self,dic,codon):
        AA = self.CodonToAmino(codon)
        thisList = dic[AA]
        freq = 0.0
        count = 0
        for entry in thisList:
            if codon == entry[0]:
                freq = entry[1]
                return (codon,count)
            count += 1

    def Translate (self):
        '''
        Takes in a list of tuple with the first value being the AA and the
        second value being its frequency in the reference genome
        '''
        finalSeq = ''
        self.sortDict(self.CBDcon)
        for entry in self.codonList:
            codon = entry[0]
            position = entry[1]
            finalSeq += self.optimize(codon,position)
        return finalSeq

    def optimize(self, codon, position):
        AA = self.CodonToAmino(codon)
        ThisList = self.CBDcon[AA]
        return ThisList[position][0]

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

'''
#Make sure the dictionaries have a .dic extension
try:
    refDic.split('.')[1]
    convDic.split('.')[1]
except IndexError:
    print("\nError: Dictionaries must end in an extension.\n")
    print(usage)
    sys.exit(2)
'''
#Create an object for the two dictionaries
mydic = read_dic(refFile,convFile)

#mydic.MakeDict(refFile,convFile)
mydic.MakeDict()


#Access a dictionary file name
mydic.sortDict(mydic.CBDcon)
mydic.sortDict(mydic.CBDref)

mydic.AnalyzeRef(seqFile)

#Access a AA within a dictionary
#protein = mydic.AnalyzeRef("ATGATCTATAAGTAA")
Translation = mydic.Translate()
OUTFILE = "Optimized-"+seqFile
FileOutput = open(OUTFILE, 'w+')
print(Translation,file = FileOutput)
''' 
Need to make sure that the sequence file exists
Add computations for the dictionary.
'''
