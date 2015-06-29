#!/usr/bin/env python3
#Name: Raymond Bryan (rvbryan)
#Group Members: Kendy Hoang (khoang3)
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
	# Second Base
	     # U             C             A             G
	#U
	'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
	'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
	'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
	'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
	#C 
	'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
	'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
	'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
	'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
	#A
	'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
	'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
	'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
	'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
	#G
	'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
	'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
	'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
	'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
                }
    def __init__ (self):
        '''In this method we will be able to analyse the sequences and then
        calculating the frequencies of each codon. '''
        #The following creates a clear table of possible codons
        self.CodonCount = {}

        for codon, AA in CodonFreq.rnaCodonTable.items():
            try:
                self.CodonCount[AA].append( {codon : 0})
            except KeyError:
                self.CodonCount[AA] = [{codon: 0}]

        #This will act as a list of the ordered frequencies
        CodonFreq = []

        #this dictionary is used to translate and transcribe the reverse
        #complement of a sequence
        self.nucs = {'A':"T", 'T':"A", 'G':"C", 'C':"G"}

        #Master list of ORFS
        self.ORFlist = []
		
        self.forwardStarts = [ [0], [0], [0] ]
		
        self.reverseStarts = [ [0], [0], [0] ]

    def readSeqRev (self, Sequence):
        '''The following function allows the reverse read to be	analysed of
        ORF's and append them to a master list.	'''
        #transcribes the reverse compliment of the orig. strand
        for nuc, newNuc in self.nucs.items():
            Comp = Sequence.replace(nuc, newNuc)

        revComp = Comp[::1]

        #Go over each nucleotide in the sequence 
        for i in range (0, len(revComp)):
            codon = revComp[i : i+3]
            #using a continuous frame will allow codon counts 
            #in between the range
            contFrame = i%3
            #Finding all start codons
            if codon in findORFs.startCodons:
                #Determines what frame the starting codon is in
                frame = i%3

                self.reverseStarts[frame].append(i)

                #save and count the codons in between once a start is filled
                #but only at the start codon frame
                if self.reverseStarts[contFrame] != []:
                    #Find what AA the codon codes for
                    Amino = CodonFreq.rnaCodonTable[codon]
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

                    #Finding the next stop and calculating the gene length
                    if codon in findORFs.stopCodons:
                        framePosition = i%3
                        #checks to see if the stop has an associated start
                        if self.reverseStarts[framePosition] != []:
                            #reset the frame        
                            self.reverseStarts[framePosition] = []

            return self.CodonCount
		
    def readSeq(self, OriginalSeq):
        '''Same as Reverse read except with the forward seq'''
        #Go over each nucleotide in the sequence 
        for i in range (0, len(OriginalSeq)):
            contFrame = i%3
            codon = OriginalSeq[i : i+3]
            #Finding all start codons
            if codon in findORFs.startCodons:
                #Determines what frame the starting codon is in
                frame = i%3
                self.forwardStarts[frame].append(i)

                #save and count the codons in between once a start is filled
                #but only at the start codon frame
                if self.reverseStarts[contFrame] != []:
                    #Find what AA the codon codes for
                    Amino = CodonFreq.rnaCodonTable[codon]
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
					
                #Finding the next stop and calculating the gene length
                if codon in findORFs.stopCodons:
                    framePosition = i%3
                    #checks to see if the stop has an associated start
                    if self.forwardStarts[framePosition] != []:
                        #reset the frame        
                        self.forwardStarts[framePosition] = []

        return (self.CodonCount)
			
    def TableMaker(self):
        '''This program  '''
	
        for AA in self.CodonCount.keys():
            #get total times the AA is coded
            for CodonANDcount in self.CodonCount[AA]:
                TotalPoss = 0 
                for Codon, Count in CodonANDcount.items():
                    TotalPoss += Count
			
        #for each codon within the AA group calc frequency
        #and replace the count for the frequency
                for i in range(0,len(self.CodonCount[AA])):
                    for Codon, Count in self.CodonCount[AA][i].items():
                        Freq = Count/TotalPoss
                        FreqPercent = float("{0:.2f}".format(Freq))
                        self.CodonCount[AA][i][Codon] = FreqPercent
			
        return (self.CodonCount)

import sequenceAnalysis
def main():
    '''Implements the Usage exception handler that can be raised from anywhere in process.'''
	 
		
    Class = CodonFreq()
	
    seqRead = sequenceAnalysis.FastAreader('tass2.fa')
 

    for header, Sequence in seqRead.readFasta():

        print (header)
        Class.readSeq(Sequence)
        Class.readSeqRev(Sequence)
        Class.TableMaker()
        print (Class.CodonCount())
		
 
if __name__ == "__main__":
	main()



		
		
	
	
	
