# April 12, 2019
# Tanmay Tanna

from __future__ import division
import sys, os, argparse, operator, numpy, pandas, fuzzysearch
from collections import Counter


## set up parser for user inputs 
parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-p', '--inFile', help='path to .fasta file.', dest='file_inFile')
required.add_argument('-o', '--outPath', help='path to output directory.', dest='path_outPath')
required.add_argument('-l', '--LBC', help='library identifier ', dest='LBC', default='null')
args = parser.parse_args()

# assign arguments to variables

inFile = str(args.file_inFile)
outPath = str(args.path_outPath)+'/'
LBC = str(args.LBC)
outName = inFile.split("/")[-1]
if outName.endswith('.fasta'):
    outName = outName[:-6]

if ('.fasta' in inFile):
    
    # open inFile for reading/writing and report file being processed
    F = open(inFile,mode='rU')
    G = open(outPath+outName+'.noLBC.fasta',mode='w')  # unique spacers based on spacer sequence only
    rawReads=0
    NonLBCReads=0

    os.system(str("echo '##################################################'"))
    os.system(str('echo '+"'"+inFile+' accepted for processing'+"'"))


   
    for L in F: # loop through reads in file
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            rawReads+=1
            continue
        L=L.strip()
            

        # Identify LBC within file
        if 'null' not in LBC:
            findLBC = L.find(LBC)
            if findLBC==-1:
                numMismatches = int(round(len(LBC)*0.1))
                LBCcoord = fuzzysearch.find_near_matches(LBC, L, max_l_dist = numMismatches)
            else:
                LBCcoord = [[findLBC, findLBC+len(LBC)]]

            if not LBCcoord:
                NonLBCReads+=1
                G.write(readName+'\n'+L+'\n')
                continue
            
        
    os.system(str('echo '+"'"+'*'+'\t'+' sampleName'+'\t'+' rawReads'+'\t'+'NonLBCReads'+"'"))
    os.system(str('echo '+"'"+'@'+'\t'+str(outPath+outName)+'\t'+str(rawReads)+'\t'+str(NonLBCReads)+"'"))



    F.close()
    G.close()
    