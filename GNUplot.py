'''
This program will create a file that can be used in GNUplot.

The first column will be the peptide number. Every column after that
will be the frequency for each codon.

Usage: 
    python3 GNUplot.py Bias_Table Nucleotide_Seq

'''

class GNUmaker:

    def __init__(self,proteinSeq,freqSeq,SSpred,Thousand):
        self.proteinSeq = proteinSeq
        self.freqSeq = freqSeq
        self.SSpred = SSpred
        self.Thousand = Thousand
        self.reverseTrans = {'C':1,'H':2,'E':4}

    def makeTable(self):
        for count,b  in enumerate(zip(self.proteinSeq,self.freqSeq)):
            AA = b[0]
            freq = b[1]
            thousand = self.Thousand[count]
            reciprocal = 1/int(freq) * thousand
            # Need to subtract 1 from count on the final seq
            try:
                SS = self.reverseTrans[self.SSpred[count][0]]
            except IndexError:
                continue
            print("%d  %s  %s  %d  %.2f  %.2f" % (count+1,AA,freq, SS, thousand, reciprocal))


def main():
    pass

if __name__ == "__main__":
    main()
