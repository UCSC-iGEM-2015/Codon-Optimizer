#!/usr/bin/env python3

# Author: David L. Bernick
# Editted by: Jairo Navarro


########################################################################
# CommandLine
########################################################################
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
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'The FOCUS program is used to view rare frequency codons using a codon usage table and ', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] CodonBiasTable NucleicAcidSeq'
                                             )
        self.parser.add_argument('CodonFile', action = 'store', help='Codon Bias Table in GCG format')
        self.parser.add_argument('NucFile', action = 'store', help='Nucleic Acid sequence')
        self.parser.add_argument('SecStruct', action='store',default=None, nargs='?',help='Secondary Structure Prediction in I-TASSER format')
        self.parser.add_argument('OutFile', action = 'store', default=None,nargs='?',help='Output file name')
#        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
#        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), action = 'store', help='minimum Gene length')
#        self.parser.add_argument('-s', '--start', action = 'append', nargs='?', help='start Codon') #allows multiple list options
#        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   

def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.  

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(myCommandLine)
###### replace the code between comments.
        # myCommandLine.args.inFile has the input file name
        # myCommandLine.args.outFile has the output file name
        # myCommandLine.args.longestGene is True if only the longest Gene is desired
        # myCommandLine.args.start is a list of start codons
        # myCommandLine.args.minGene is the minimum Gene length to include
        #
    print (myCommandLine.args) # print the parsed argument string .. as there is nothing better to do

    if True:
        print ('First Arg is', str(myCommandLine.args.CodonFile) )
        print (type(myCommandLine.args.CodonFile))
    else :
        pass
#######
    
if __name__ == "__main__":
    main()


