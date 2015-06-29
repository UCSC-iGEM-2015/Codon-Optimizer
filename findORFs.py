#!/usr/bin/env python3
#Name: Raymond Bryan (rvbryan)
#Group Members: Kendy Hoang (khoang3)
'''
The following program will take and open the file it has been given 
in its directory. The program will then analyse and differentiate each
sequence header to its corresponding nucleotide sequences. It will
then take the sequence and go through the sequence taking all possible genes out 
that are over 100 nucleotides in length.


user inputs: File named "tass2.fa.fa" in same directory
program outputs: Each gene above the 100 nucleotides according to frame, position,
and gene length. In addition to printing this data it will continue to write the data
to a txt file and saved in the same directory as 'tass2ORFdata-ATG-100.txt'

'''

###############Pseudo Code: Start################################
'''
1. Taking the file given and extracting each corresponding header/seq.(FastAreader)

2. Calculate the reverse compliment of the seq given. 

3. Take the reverse complement of the seq and locate all possible 
        start codons
        
4. Determine Frame

5. Locate a stop codon in connection with the previous start and frame.

6. Do the same with forward strand

7. Set boundary conditions:

'''
###############Pseudo Code:END###################################


######################################################################
 #CommandLine
######################################################################

import sequenceAnalysis #module import for fasta reader
imprt
                
class findORFs:
        '''FindOrfs is our class that contains the functions that mediate
        orf finding. 
        '''
        startCodons = ['ATG', 'CTG', 'GTG']
        stopCodons = ['TGA','TAG','TAA']

        def __init__ (self):
                '''The init method intiates the dictionaries and lists we will use for the entire program '''
                
                #this dictionary is used to translate and transcribe the reverse
                #complement of a sequence
                self.nucs = {'A':"T", 'T':"A", 'G':"C", 'C':"G"}

                #Master list of ORFS
                self.ORFlist = []
                
                self.forwardStarts = [ [0], [0], [0] ]
                
                self.reverseStarts = [ [0], [0], [0] ]
        
        def readSeqRev (self, Sequence):
                '''The following function allows the reverse read to be
                analysed of ORF's and append them to a master list. 
                '''
                #transcribes the reverse compliment of the orig. strand
                for nuc, newNuc in self.nucs.items():
                        Comp = Sequence.replace(nuc, newNuc)
                
                revComp = Comp[::1]
                
                #Go over each nucleotide in the sequence 
                for i in range (0, len(revComp)):
                        codon = revComp[i : i+3]
                        #Finding all start codons
                        if codon in findORFs.startCodons:
                                #Determines what frame the starting codon is in
                                frame = i%3
                                self.reverseStarts[frame].append(i)

                        #Finding the next stop and calculating the gene length
                        if codon in findORFs.stopCodons:
                                framePosition = i%3
                                #checks to see if the stop has an associated start
                                if self.reverseStarts[framePosition] != []:
                                
                                        startPositon = len(revComp) - (i + 2)
                                        stopPosition = len(revComp) - self.reverseStarts[framePosition][0]
                                        geneLength = stopPosition - startPositon
                                        
                                        if geneLength > 100:
                                                #add to the list of ORFs
                                                self.ORFlist.append([startPositon +1 , stopPosition + 1, geneLength + 1, -(framePosition + 1)])
                                        #reset the frame        
                                        self.reverseStarts[framePosition] = []
                                        
                return (self.ORFlist)
                
        def readSeq(self, OriginalSeq):
                '''Same as Reverse read except with the forward seq'''
                #Go over each nucleotide in the sequence 
                for i in range (0, len(OriginalSeq)):
                        codon = OriginalSeq[i : i+3]
                        #Finding all start codons
                        if codon in findORFs.startCodons:
                                #Determines what frame the starting codon is in
                                frame = i%3
                                self.forwardStarts[frame].append(i)

                        #Finding the next stop and calculating the gene length
                        if codon in findORFs.stopCodons:
                                framePosition = i%3
                                #checks to see if the stop has an associated start
                                if self.forwardStarts[framePosition] != []:
                                        
                                        startPositon = self.forwardStarts[framePosition][0]
                                        stopPosition = i + 2
                                        geneLength = stopPosition - startPositon
                                        
                                        if geneLength > 100:
                                                #add to the list of ORFs#add 1
                                                self.ORFlist.append([startPositon +1 , stopPosition +1, geneLength + 1, (framePosition + 1)])
                                        #reset the frame        
                                        self.forwardStarts[framePosition] = []
                                        
                return (self.ORFlist)


class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        import argparse
        
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does',
                                              epilog = 'Program epilog - some other stuff you feel compelled to say',
                                              add_help = True, #default is True
                                              prefix_chars = '-',
                                              usage = '%(prog)s [options] -option1[default] <input >output')
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=True, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= range(0, 2000), action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
                self.args = self.parser.parse_args()
        else :
                self.args = self.parser.parse_args(inOpts)

######################################################################
#Main
#Here is the main program
 

######################################################################
   

def main():
        '''Implements the Usage exception handler that can be raised from anywhere in process.'''
         
        orfOutput = open('./tass2ORFdata-ATG-100.txt', "w")
        myCommandLine = CommandLine([ 'tass2.fa',
                                      'tass2ORFdata-ATG-100.txt',
                                      '--longestGene',
                                      '--start=ATG',
                                      '--minGene=100'])


                
        ORF = findORFs() 
        seqRead = sequenceAnalysis.FastAreader('tass2.fa.fa')


        for header, Sequence in seqRead.readFasta():

                print (header)
                print(header, file=orfOutput)
                ORF.readSeq(Sequence)
                ORF.readSeqRev(Sequence)
                for i in range(len(ORF.ORFlist)):
                        ORF.ORFlist.sort(key = lambda entry: entry[2], reverse = True)   #sort each orf by decreasing length

                #print: frame, start position, end position, length
                        print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(ORF.ORFlist[i][3], ORF.ORFlist[i][0], ORF.ORFlist[i][1], ORF.ORFlist[i][2]))
                        print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(ORF.ORFlist[i][3], ORF.ORFlist[i][0], ORF.ORFlist[i][1], ORF.ORFlist[i][2]), file = orfOutput)
 
if __name__ == "__main__":
    main()

