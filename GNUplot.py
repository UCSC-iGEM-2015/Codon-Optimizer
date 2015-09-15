'''
Author: Jairo Navarro
Date: July 23, 2015

Description:
This program will create a file that can be used in a 
GNUplot script. Currently, the program only contains the class that is
used in the FOCUS python program, and does not work on its own.

Formatting:
The first column will be the amino acid position.
    second column will be the amino acid
    third column will be the number 1, 2 or 4
        -Coil: 1
        -Helix: 2
        -Sheet: 4
    fourth column will be the percent per thousand for the codon

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
            
            try:
                SS = self.reverseTrans[self.SSpred[0][count]]
            except IndexError:
                continue
            print("%d  %s  %s  %d  %.2f" % (count+1,AA,freq, SS, thousand))


def main():
    pass

if __name__ == "__main__":
    main()
