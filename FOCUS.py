'''
Name: Jairo Navarro
Date: July 19, 2015

The following program is a command line program for the use of codon optimization. Currently, the program only views where 'rare' codons are likely to be in the nucleic acid sequence, but has noise that needs to be removed. So take the codon 'score' with a grain of salt. 

This program takes in an organism's codon bias and displays a numerical
representation of the relative codon frequency for each amino acid 
residue in a particular sequence from that organism. It also takes in an
I-TASSER format secondary structure prediction of the same sequence ran through PSSpred.


Usage: 

> python3 FOCUS.py CodonBiasTable Sequence SSpred

Where: 
FOCUS.py is the program name
CodonBiasTable is the GCG codon bias table for the reference organism.
Sequence is the nucleic acid sequence for the desired protein
SSpred is the secondary structure prediction file in I-TASSER format. 

Output:
AA sequence with relative codon frequencies, SS prediction, and confidence 
displayed below in a text file

'''
class FOCUS:

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
        # Convert the codon to an Amino Acid
        return self.dnaCodonTable[codon]

    def MakeDict (self):
        # Open the conversion codon bias tables
        gcgFile = self.OpenFile(self.GCG)
        # Read line by line
        for line in gcgFile:
            # Ingore newlines and whitespace between data in file
            if line.strip() == '':
                continue
            # Extract information from the codon bias table
            codon = line.strip().split()[1]
            AA = self.CodonToAmino(codon)
            freq = float(line.split()[4])
            thousand = float(line.split()[3])
            # For the current amino acid, append the list to the dictionary
            try:
                self.CodonMap[AA].append([codon,freq,thousand])
            except:  
                #If it is the firt codon for the amino acid, add the key and value
                self.CodonMap[AA] = [[codon,freq,thousand]]
        #Sort the codon dictionary if necessary
        self.CodonMap = self.sortDict(self.CodonMap)

    def sortDict(self, dic):
        # Sort codons per AA by frequency
        # Highest frequency first

        from operator import itemgetter
        for entry,value in dic.items():
            sortedList = sorted(value, key=itemgetter(1), reverse=True)
            dic[entry] = sortedList
        return dic

    def getScore(self,codon):
        '''Returns the relative codon score for a paticular amino acid'''

        AA = self.CodonToAmino(codon)
        codonList = self.CodonMap[AA]
        position = 0
        # Go through the list until codon is found
        for entry in codonList:
            if codon == entry[0]:
                # Break once position is found
                break
            position += 1
        return position + 1

    def printer(self,prot,values):
        '''Format the protein with the codon frequency'''
        #Counter for the entire protein sequence length
        count = len(prot)
        #Start at position 0
        start = 0
        #Have a step of 71, that is print 71 characters per line
        step = 71
        #Final string container
        finalSeq = ''

        #While the counter is greater than 0, add to the format string
        while(count>0):
            #Add the protein sequence
            finalSeq += "%s\n" % prot[start:start+step]
            #Add the codon score
            finalSeq += "%s\n\n" % values[start:start+step]

            #Decrease the counter by the step
            count -= step
            #Increase the starting position by the step
            start += step

        #Return the final formatted string
        return finalSeq


    def SSprinter(self,prot,values,SSlist):
        '''
        This code is for printing when given the secondary structure
        file. Similar to the above function but with added features.
        '''

        #Counter for the entire protein sequence length
        count = len(prot)
        #Start at position 0
        start = 0
        #Have a step of 71, that is print 71 characters per line
        step = 71
        #Get the secondary structure string from the SSlist
        SecondaryStructure = SSlist[0]
        #Get the confidence for the Secondary Structure from SSlist
        Confidence = SSlist[1]
        #Final string container
        finalSeq = ''

        #While the counter is greater than 0, add to the format string
        while(count>0):
            #Add the protein sequence
            finalSeq += "%s\n" % prot[start:start+step]
            #Add the codon score
            finalSeq += "%s\n" % values[start:start+step]
            #Add the secondary structure
            finalSeq += "%s\n" % SecondaryStructure[start:start+step]
            #Add the confidence of Secondary Structure and extra space
            finalSeq += "%s\n\n" % Confidence[start:start+step]
            #Decrease the counter by the step
            count -= step
            #Increase the starting position by the step
            start += step

        #Return the final formatted string
        return finalSeq


'''
Main program to be ran
'''
def main():
    '''Start the program with some imports'''

    # fastaReader to read FASTA files
    import fastaReader as fRead
    # parseTasser to parse the Secondary structure file
    import parseTasser as pT
    # commandLine to parse command line options
    import commandLine

    #Command line object with given options
    myCommandLine = commandLine.CommandLine()

    #Codon Bias Table for the reference organism
    REF_ORG = myCommandLine.args.CodonFile
    #Nucleotide Sequence for desired protein
    FASTA = myCommandLine.args.NucFile
    #Secondary Structure for the desired protein
    Align = myCommandLine.args.SecStruct
    #Output file name if specified
    OutFile = myCommandLine.args.OutFile
    #If no output file name was specified, just add focus to the 
    # begining of the FASTA file
    if OutFile == None:
        OutFile = "FOCUS-" + FASTA

    # File that the program writes to...
    FileOutput = open(OutFile,'w+')

    #Parse secondary structure file if given
    if Align != None:
        thisParser = pT.parseTASSER(Align)
        #SSlist is a list with two entries
        # SSlist[0] = string for the entire secondary structure
        # SSlist[1] = string for the confidence of secondary structure
        SS_list = thisParser.returnList()

    # containers for FASTA reader and Optimizer
    thisReader = fRead.FastAreader(FASTA)
    Recoder = FOCUS(REF_ORG)

    # Make and sort codon bias tables by frequency
    Recoder.MakeDict()

    # Read through the nucleotide FASTA File
    for head,seq in thisReader.readFasta():
        '''String containers'''
        # Score for the given codon
        value = ''
        # Protein sequence
        proteinSeq = ''
        #Percent per thousand list
        Thousand = []

        # For each codon in the FASTA file
        for nuc in range(0, len(seq), 3):
            # Save the codon
            codon = seq[nuc:nuc+3]
            # Get the Amino acid for the given codon
            AA = Recoder.CodonToAmino(codon)
            # Add the amino acid to the protein sequence
            proteinSeq += AA
            #FOCUS codon score
            codonScore = Recoder.getScore(codon)
            # Get the frequency per thousand from the codon bias table dictionary
            FrequencyPerThousand = Recoder.CodonMap[AA][codonScore-1][2]
            # Append the the thousand list
            Thousand.append(FrequencyPerThousand)
            # Add the relative codon freq to the score string
            value += str(codonScore)

        # Print the header to the file
        print('>' + head, file = FileOutput)
        # Check to see if secodnary structure file was given
        if Align != None:
            # If true, use the Secondary structure 'printer'
            FormattedSeq = Recoder.SSprinter(proteinSeq,value,SS_list)
            
            #import the GNUplot python program
            #This is used to print the data in a format that can be
            # used in a GNUplot script. Can be removed if not needed.
            # This does not write to the file, but printed to stdout
            import GNUplot as GP
            temp = GP.GNUmaker(proteinSeq,value,SS_list,Thousand)
            temp.makeTable()
        else:
            # Otherwise, use the default 'printer'
            FormattedSeq = Recoder.printer(proteinSeq,value)


        # Write the formatted sequence from the 'printer' to the file
        print(FormattedSeq, file = FileOutput)

        

if __name__ == "__main__":
    main()

