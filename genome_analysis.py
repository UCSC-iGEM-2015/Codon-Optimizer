import sys
import read_genome as read

fileName = sys.argv[0]
fileH = sys.argv[1]

print("Program Name: %s" % fileName)
print("FASTA File: %s" % fileH)
newname = fileH.split(".")
newname = newname[0] + ".dic"
print(newname)
thisSeq = ''
mygenome = read.read_genome(fileH)
count = 0
for head,seq in mygenome.ReadFile():
    print(head)
    if(count>5): 
        break
    thisSeq += seq
