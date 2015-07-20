'''
Name: Jairo Navarro, John Poncini
Date: July 19, 2015

This program takes in an organism's codon bias and displays a numerical
representation of the relative codon frequency for each amino acid 
residue in a particular sequence from that organism. It also takes in an
I-TASSER format SS prediction of the same sequence ran through PSSpred.

The following is the secondary structure prediction in I-TASSER format
where the 3th column is the SS prediction (1: coil; 2: helix; 4: strand)
and the 4th column is confidence of the prediction (0: low; 9: high).

Usage: 

> python3 BiasAlign.py CodonBiasTable Sequence SSpred
> Enter codon bias table in GCG format:         codonbias.txt
> Enter genomic sequence in FASTA format:       seq.fa
> Enter PSSpred prediction in I-TASSER format:  seq.dat

Output:
AA sequence with relative codon frequencies, SS prediction, and confidence 
displayed below in a text file

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

    def __init__ (self, GCG):
        '''
        Containers for the reference dictionary
                           converstion dictionary
        '''
        # Save the file names
        self.GCG = GCG
        #CodonBiasDicionary
        self.CodonMap = {}


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

    def MakeDict (self):
        # Open the conversion codon bias tables
        gcgFile = self.OpenFile(self.GCG)
        # Read line by line
        for line in gcgFile:
            # Codon from the Codon Bias Table
            codon = line.split()[1]
            AA = self.CodonToAmino(codon)
            freq = float(line.split()[4])
            
            try:
                self.CodonMap[AA].append((codon,freq))
            except:  
                self.CodonMap[AA] = [(codon,freq)]
        self.CodonMap = self.sortDict(self.CodonMap)
        print(self.CodonMap)
    def sortDict(self, dic):
        # Sort codons per AA by frequency
        # Highest frequency first

        from operator import itemgetter
        for entry,value in dic.items():
            sortedList = sorted(value, key=itemgetter(1), reverse=True)
            dic[entry] = sortedList
        return dic

    def recode(self,codon):
        AA = self.CodonToAmino(codon)
        codonList = self.CodonMap[AA]
        print(codonList)
        position = 0
        # Go through the list until codon is found
        for entry in codonList:
            if codon == entry[0]:
                # Break once position is found
                break
            position += 1
        return position + 1

    def printer(self,prot,values):
        outerCount = 0
        innerCount = 0
        tracker = 0
        finalSeq = ''
        for AA in prot:
            if outerCount >= 71:
                outerCount = 0 
                finalSeq += "\n"
                for num in range(tracker,len(values)):
                    if innerCount >= 71:
                        innerCount = 0
                        finalSeq += "\n"
                        break
                    finalSeq += "%s" % values[tracker]
                    innerCount += 1
                    tracker += 1
                finalSeq += "\n"
            finalSeq += "%s" % AA
            outerCount += 1
        finalSeq += "\n"
        for num in range(tracker, len(values)):
            finalSeq += "%s" % values[tracker]
            tracker += 1
        return finalSeq

def main():
    import fastaReader as fRead
    import sys
    #####################################################
    '''Replace with CommandLine class'''

    #First file given in command line after program name
    REF_ORG = sys.argv[1]
    # FASTA file
    FASTA = sys.argv[2]
    OutFile = "Aligned-" + FASTA
    
    #####################################################

    FileOutput = open(OutFile,'w+')

    thisReader = fRead.FastAreader(FASTA)
    Recoder = optimizer(REF_ORG)

    # Make and sort codon bias tables by frequency
    Recoder.MakeDict()

    for head,seq in thisReader.readFasta():
        value = ''
        proteinSeq = ''

        for nuc in range(0, len(seq), 3):
            codon = seq[nuc:nuc+3]
            AA = Recoder.CodonToAmino(codon)
            proteinSeq += AA
            value += str(Recoder.recode(codon))
        
        print('>' + head, file = FileOutput)
        FormattedSeq = Recoder.printer(proteinSeq,value)
        '''Print on own version formatted'''
        #print(proteinSeq)
        #print(value)
        '''Print on own version formatted'''
        
        print(FormattedSeq, file = FileOutput)

    print("\nSequence saved as: \n\t*%s\n" % OutFile)
        

if __name__ == "__main__":
    main()
