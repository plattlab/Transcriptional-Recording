from __future__ import division
import sys, os, argparse, operator, numpy



## set up parser for user inputs 
parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-i', '--inFile', help='path to .sam file.', dest='file_inFile')
required.add_argument('-o', '--outFile', help='path to output file', dest='path_outFile')


args = parser.parse_args()

# assign arguments to variables

inFile = str(args.file_inFile)
outFile = str(args.path_outFile)

# identify and process files with the terms below
if ('.sam' in inFile):
    
    # open inFile for reading/writing and report file being processed
    F = open(inFile, mode='rU')
    G = open(outFile, mode='w')  
    alignments = {}

   
    for L in F: # loop through reads in file
        if '@' in L: # defline, skip for processing but save read name
            G.write(L)
            continue
        start_coord = L.split('\t')[3]
        seq = L.split('\t')[9]

        if start_coord in alignments:
            if len(seq) in alignments[start_coord]:
                continue
            else:
                alignments[start_coord].append(len(seq))
                G.write(L)
        else:
            alignments[start_coord]=[len(seq)]
            G.write(L)

    F.close()
    G.close()

else:
    print('invalid input file!')
