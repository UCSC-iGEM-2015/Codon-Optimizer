'''
Python program to parse I-TASSER secondary structure format files
The program as currently written returns the 3rd and 4th column
'''

class parseTASSER:

    def __init__(self,TASSER):

        #Saves file name
        self.TASSER = TASSER
        self.trans = {1:'C',2:'H',4:'E'}

    def OpenFile (self,thisFile):
        try:
            return open(thisFile,'r')
        except FileNotFoundError:
            print("\n%s does not exist.\n" % thisFile)

    def returnList (self):
        returnList = []
        #Open the secondary structure file
        tassFile = self.OpenFile(self.TASSER)
        # Read line by down
        for line in tassFile:
            line = line.split()
            thisTUPLE = (self.trans[int(line[2])],int(line[3]))
            # thisTUPLE = (Letter,Confidence)
            returnList.append(thisTUPLE)
        return returnList

def main():
    import sys

    FILE = sys.argv[1]
    parser = parseTASSER(FILE)
    LIST = parser.returnList()
    for item in LIST:
        print(item)

if __name__ == "__main__":
    main()

