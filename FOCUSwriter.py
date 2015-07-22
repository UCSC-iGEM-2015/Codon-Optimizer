'''
This program takes in a nucleotide sequence and prints the 
amino acid sequence above the codons. The codons are also 
separated by a space. 
'''

class writer:

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
    dnaCodonTable = {key.replace('U','T'): value for key, value in rnaCodonTable.items()}
    
    def __init__ (self, FOCUS, biasFile):
        self.FOCUS = FOCUS
        self.biasFile = biasFile
        self.CodonMap = {}

    def CodonToAmino(self, codon):
        # Change DNA to RNA and return and AA
        return self.dnaCodonTable[codon]

    def OpenFile(self, seqFile):
        '''Open the sequence file and ensure that it exists.'''
        try:
           return open(seqFile,'r')
        except FileNotFoundError:
           print("\n%s does not exist.\n" % seqFile)
           sys.exit(3)

    def MakeDict (self):
        # Open the conversion codon bias tables
        gcgFile = self.OpenFile(self.biasFile)
        # Read line by line
        for line in gcgFile:
            if (line.startswith('\n')):
                continue
            # Codon from the Codon Bias Table
            codon = line.split()[1]
            AA = self.CodonToAmino(codon)
            freq = float(line.split()[4])
            thousand = float(line.split()[3])
            try:
                self.CodonMap[AA].append([codon,freq,thousand])
            except:  
                self.CodonMap[AA] = [[codon,freq,thousand]]
        self.CodonMap = self.sortDict(self.CodonMap)

    def sortDict(self, dic):
        # Sort codons per AA by frequency
        # Highest frequency first

        from operator import itemgetter
        for entry,value in dic.items():
            sortedList = sorted(value, key=itemgetter(1), reverse=True)
            dic[entry] = sortedList
        return dic

    def readFile(self):
        thisFile = self.OpenFile(self.FOCUS)
        header = ''
        proteinSeq = ''
        value = ''
        count = 0
        for line in thisFile:
            print(line)
            
            if (line.startswith('>')):
                header = line.split('\n')[0]
                continue
            elif (line.startswith('\n')):
                continue
            else:
                if count %2 == 0:
                    proteinSeq += line.split('\n')[0]
                if count %2 == 1:
                    value += line.split('\n')[0]
                count += 1
        thisList = [header,proteinSeq,value]
        return thisList

    def rewrite(self,protein,freq):
        nucSeq = ''
        count = 0
        for AA, freq in zip(protein,freq):
            if count >= 23:
                nucSeq += '\n'
                count = 0
            codon = self.CodonMap[AA][int(freq)-1][0]    
            nucSeq += codon
            count += 1
        return nucSeq

def main():
    # Some imports
    import fastaReader as fRead
    import sys

    #First file given in command line after program name
    BIAS = sys.argv[1]
    # FASTA file
    FASTA = sys.argv[2]
    outfile = "Optimized-%s" % FASTA

    Recoder = writer(FASTA,BIAS)
    Recoder.MakeDict()
    DATA_LIST = Recoder.readFile()
    header = DATA_LIST[0]
    proteinSeq = DATA_LIST[1]
    values = DATA_LIST[2]
    FILEoutput = open(outfile,'w+')

    print(header,file=FILEoutput)
    nucleotideSeq = Recoder.rewrite(proteinSeq,values)
    print(nucleotideSeq,file=FILEoutput)

if __name__ == "__main__":
    main()
