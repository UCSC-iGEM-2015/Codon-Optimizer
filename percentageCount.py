#!/usr/bin/env python3
#Name: Raymond Bryan (rvbryan)
#Date: 6/30/2015
'''
The following program will ask the user for
a protein sequence. The program will then change whatever
input given to a pure AA string. Once in this string, the
program will exclude any anomaly not one of the 20 Aminoacids.
user inputs: Protein sequence
program outputs: Number of Amino Acids, Molecular Weight,
molar Extinction coefficient, mass Extinction coefficient,
Theoretical pI, and Amino acid composition by percentage.
'''


class CodonFreq:
    """
    This class will give the frequentcies of each condon in the form a 
    dictionary. The dictionary will be coded from a check system that is provided 
    """

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
        '''In this method we will be able to analyse the sequences and then
        calculating the frequencies of each codon. '''
        #The following creates a clear table of possible codons
        self.CodonCount = {}

        for codon, AA in CodonFreq.dnaCodonTable.items():
            try:
                self.CodonCount[AA].append( {codon : 0})
            except KeyError:
                self.CodonCount[AA] = [{codon: 0}]

        
        self.CodonFrequen = {}
        for codon, AA in CodonFreq.dnaCodonTable.items():
            try:
                self.CodonFrequen[AA].append( {codon : 0})
            except KeyError:
                self.CodonFrequen[AA] = [{codon: 0}]
                
        #this dictionary is used to translate and transcribe the reverse
        #complement of a sequence
        self.nucs = {'A':"T", 'T':"A", 'G':"C", 'C':"G"}

        #Master list of ORFS
        self.ORFlist = []
		
        self.forwardStarts = [ [0], [0], [0] ]
		
        self.reverseStarts = [ [0], [0], [0] ]

        self.startCodons = ['ATG', 'CTG', 'GTG']
        self.stopCodons = ['TGA','TAG','TAA']

    def readSeqRev (self, Sequence):
        '''The following function allows the reverse read to be	analysed of
        ORF's and append them to a master list.	'''
        #transcribes the reverse compliment of the orig. strand
        for nuc, newNuc in self.nucs.items():
            Comp = Sequence.replace(nuc, newNuc)

        revComp = Comp[::1]

        #Go over each nucleotide in the sequence 
        for i in range (0, len(revComp),3):
            codon = revComp[i : i+3]
            
            #Find what AA the codon codes for
            try:
                Amino = CodonFreq.dnaCodonTable[codon]
            except KeyError:
                continue
            #Find how many codons code for that AA
            possibleCodons = len(self.CodonCount[Amino])
            
            #search possible Codons for the existing count
            for position in range(0,possibleCodons):
                #increment count of that codon by 1
                try:
                    exCount = self.CodonCount[Amino][position][codon]
                    exCount += 1
                    self.CodonCount[Amino][position][codon] = exCount
                except KeyError:
                    #try the next one
                    continue
        
        return self.CodonCount
		
    def readSeq(self, OriginalSeq):
        '''Same as Reverse read except with the forward seq'''
        #Go over each nucleotide in the sequence 
        for i in range (0, len(OriginalSeq), 3):
            codon = OriginalSeq[i : i+3]

            #Find what AA the codon codes for
            try:
                Amino = CodonFreq.dnaCodonTable[codon]
            except KeyError:
                continue
            #Find how many codons code for that AA
            possibleCodons = len(self.CodonCount[Amino])

            #search possible Codons for the existing count
            for position in range(0,possibleCodons):
                #increment count of that codon by 1
                try:
                    exCount = self.CodonCount[Amino][position][codon]
                    exCount += 1
                    #print (exCount)
                    self.CodonCount[Amino][position][codon] = exCount
                    #print ('total'+str(self.CodonCount[Amino][position][codon]))
                except KeyError:
                    #try the next one
                    continue
					
        return (self.CodonCount)
			
    def TableMaker(self):
        '''This program  '''
        Total = []
        AAcount = {}
        i=0
        for AA in self.CodonCount.keys():

            #get total times the AA is coded
            for i in range(0, len(self.CodonCount[AA])):
                for codon in self.CodonCount[AA][i]:
                    Total.append(self.CodonCount[AA][i][codon])
                TotalPoss = sum(Total)
                AAcount[AA] = TotalPoss
            Total = []
            #print (AAcount)
        #for each codon within the AA group calc frequency
        #and replace the count for the frequency
        for AA in self.CodonCount.keys():
            for i in range(0,len(self.CodonCount[AA])):
                for codon in self.CodonCount[AA][i]:
                    number = self.CodonCount[AA][i][codon]
                    
                    try:
                        Freq = number/(AAcount[AA])
                        #print (Freq)
                        FreqPercent = float("{0:.2f}".format(Freq))
                        self.CodonFrequen[AA][i][codon] = FreqPercent
                    except ZeroDivisionError:
                        FreqPercent = 0


                
			
        return (self.CodonFrequen)

import sequenceAnalysis
def main():
    '''Implements the Usage exception handler that can be raised from anywhere in process.'''
    import sys
    usage = "\nUsage: \npython3 %s FASTA.fa\n" % sys.argv[0]
    if len(sys.argv) != 2:
        print("Please enter a single FASTA file.")
        print(usage)
        sys.exit(3)
    OUTFILE = sys.argv[1].split(".")[0] + ".dic"
    FASTA = sys.argv[1]
    Class = CodonFreq()
    FileOutput = open(OUTFILE, 'w+')
    seqRead = sequenceAnalysis.FastAreader(FASTA)
    print('wait for it......')
 

    for header, Sequence in seqRead.readFasta():

        #print (header)
        Class.readSeq(Sequence)
        Class.readSeqRev(Sequence)
        Class.TableMaker()

    print('wait for it......')

    for AA in Class.CodonFrequen.keys():
            for i in range(0,len(Class.CodonFrequen[AA])):
                for codon in Class.CodonFrequen[AA][i]:
                    #print: codon, AA, Frequency, number

                    print('{:s} {:s} {:.2f} {:.0f}'.format(codon, AA, Class.CodonFrequen[AA][i][codon], Class.CodonCount[AA][i][codon]))
                    print('{:s} {:s} {:.2f} {:.0f}'.format(codon, AA, Class.CodonFrequen[AA][i][codon], Class.CodonCount[AA][i][codon]), file = FileOutput)
    #print (Class.CodonCount)
    #print (Class.CodonFrequen)
		
 
if __name__ == "__main__":
	main()
