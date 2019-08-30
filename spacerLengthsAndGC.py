# August 15, 2019
# Tanmay Tanna
#k="$(ls *unique* | tr "\n" " ")"; python3 ../../../scripts/spacerLengthsAndGC.py -i $k -o spacerStats -n $i

from __future__ import division
import sys, os, argparse, operator, numpy, pandas, fuzzysearch
from collections import Counter


## set up parser for user inputs 
parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-i','--inFiles',  nargs='+', help='array of fasta files with path', required=True)
required.add_argument('-o', '--outPath', help='path to output directory.', dest='path_outPath')
required.add_argument('-n', '--outName', help='name of output file', dest='outName')
args = parser.parse_args()

# assign arguments to variables

inFileList= args.inFiles
outPath = str(args.path_outPath)+'/'
os.makedirs(outPath, exist_ok=True)
outName = str(args.outName)
G = open(outPath+outName+'.spacerLengths.txt',mode='w')
H = open(outPath+outName+'.spacerGC.txt',mode='w')
spacerLengths = {}
spacerGC = {}

for inFile in inFileList:
    if ('.fasta' in inFile):
        
        # open inFile for reading/writing and report file being processed
        F = open(inFile,mode='rU')
        
        os.system(str("echo '##################################################'"))
        os.system(str('echo '+"'"+inFile+' accepted for processing'+"'"))
       
        for L in F: # loop through reads in file
            if '>' in L: # defline, skip for processing but save read name
                continue
            L=L.strip()
            count = Counter(L)
            length=len(L)
            if length in spacerLengths:
                spacerLengths[length]+=1
            else:
                spacerLengths[length]=1
            gc = (count['G']+count['C'])*100/len(L)
            if gc in spacerGC:
                spacerGC[gc]+=1
            else:
                spacerGC[gc]=1

        F.close()

for length in spacerLengths:
    G.write(str(length)+'\t'+str(spacerLengths[length])+'\n')

for gc in spacerGC:
    H.write(str(gc)+'\t'+str(spacerGC[gc])+'\n')
         
G.close()
H.close()
    