# April 12, 2019
# Tanmay Tanna

from __future__ import division
import sys, os, argparse, operator, numpy, pandas, fuzzysearch
from collections import Counter


## set up parser for user inputs 
parser = argparse.ArgumentParser()

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-p', '--inFile', help='path to .fasta file.', dest='file_inFile')
required.add_argument('-o', '--outPath', help='path to output directory.', dest='path_outPath')
required.add_argument('-d1', '--drOne', help='first partial CRISPR repeat ', dest='dr_sequence1')
required.add_argument('-d2', '--drTwo', help='second partial CRISPR repeat ', dest='dr_sequence2')
required.add_argument('-df', '--drFull', help='second partial CRISPR repeat ', dest='dr_sequenceFull')


## optional (defaults provided)

optional.add_argument('-l', '--LBC', help='library identifier ', dest='LBC', default='null')
optional.add_argument('-s', '--outName', help='path to output directory.', dest='file_outName', default='null')
optional.add_argument('-a', help='Number of allowed mismatches in the first partial CRISPR repeat. Default=2', type=int, dest='m1', default=2)
optional.add_argument('-b', help='Number of allowed mismatches in the second partial CRISPR repeat. Default=3', type=int, dest='m2', default=3)
optional.add_argument('-m', help='Minimum spacer size. Default=30', type=int, dest='min', default=25)
optional.add_argument('-n', help='Maximum spacer size. Default=55', type=int, dest='max', default=56)
optional.add_argument('-sMin', help='Minimum stagger length according to primer design. Default=0', type=int, dest='min_stagger', default=0)
optional.add_argument('-sMax', help='Maximum stagger length according to primer design. Default=8', type=int, dest='max_stagger', default=8)
optional.add_argument('--infoFile', help='boolean to generate files with GC content and length info for spacers', dest='infoFile', action='store_true')
optional.add_argument('--no-infoFile', help='boolean to generate files with GC content and length info for spacers', dest='infoFile', action='store_false')
optional.set_defaults(infoFile=False)

parser._action_groups.append(optional) 
args = parser.parse_args()

# assign arguments to variables

DR1Mismatch = int(args.m1)
DR2Mismatch = int(args.m2)
minSpacer = int(args.min)
maxSpacer = int(args.max)
inFile = str(args.file_inFile)
outPath = str(args.path_outPath)+'/'
outName = str(args.file_outName)
firstRepeat = str(args.dr_sequence1)
secondRepeat = str(args.dr_sequence2)
fullRepeat = str(args.dr_sequenceFull)
LBC = str(args.LBC)
minStagger = int(args.min_stagger)
maxStagger = int(args.max_stagger)
infoFile=args.infoFile
if outName is 'null':
    outName = inFile.split("/")[-1]
if outName.endswith('.fasta'):
    outName = outName[:-6]



 
# Function that takes in a part of the read and gives back a spacer

def editSpacer(read,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook):
    s=''
    # Find the first repeat using string search, if not found, use fuzzy string matching

    s = read[firstExpect:(firstExpect+firstRangeToLook+len(firstRepeat))]
    findDR1 = s.find(firstRepeat)
    if findDR1 ==-1:
        firstMatch = fuzzysearch.find_near_matches(firstRepeat, s, max_l_dist = DR1Mismatch)
        firstMatch = sorted(firstMatch, key=lambda x: x[2])
    else:
        firstMatch = [[findDR1, findDR1 + len(firstRepeat)]]


    if not firstMatch:   return ''   # Too many mismatches. Return empty string.

    s=''
    # Find the second repeat using fuzzy string matching 
   
    s = read[secondExpect:(secondExpect+secondRangeToLook+len(secondRepeat))]
    if len(s)>len(secondRepeat):
        findDR2 = s.find(secondRepeat)
        if findDR2 ==-1:
            secondMatch = fuzzysearch.find_near_matches(secondRepeat, s, max_l_dist = DR2Mismatch)
            secondMatch = sorted(secondMatch, key=lambda x: x[2])
        else: 
            secondMatch = [[findDR2, findDR2 + len(secondRepeat)]]
        if not secondMatch: return ''
    else:
        return ''
    
       # Too many mismatches or read too short. Return empty string.

    # If both repeats seem to have been found, return spacer
    spacerStart = firstMatch[0][1]
    spacerEnd = secondMatch[0][0] + secondExpect
    # if spacer is too short, looking for other matches for second DR in the read
    if len(secondMatch) > 1 & spacerEnd - spacerStart < minSpacer:
        if findDR1 ==-1:
            i=1
            while i < len(secondMatch) and spacerEnd - spacerStart < minSpacer:
                spacerEnd = secondMatch[i].start + secondExpect
                i+=1
        else:
            secondMatch = fuzzysearch.find_near_matches(secondRepeat, s, max_l_dist = DR2Mismatch)
            secondMatch = sorted(secondMatch, key=lambda x: x[2])
            i=0
            while i < len(secondMatch) and spacerEnd - spacerStart < minSpacer:
                spacerEnd = secondMatch[i].start + secondExpect
                i+=1

    spacer = read[spacerStart:spacerEnd]

    # no spacer if out of bounds
    if len(spacer) > maxSpacer: return ''
    if len(spacer) < minSpacer: return ''
    
    return spacer

# identify and process files with the terms below
if ('.fasta' in inFile):
    
    # open inFile for reading/writing and report file being processed
    F = open(inFile,mode='rU')
    G = open(outPath+outName+'.unique.fasta',mode='w')  # unique spacers based on spacer sequence only
    I = open(outPath+outName+'.doubleAcquisitions.fasta',mode='w') 
    J = open(outPath+outName+'.doubleAcquisitions.paired.fasta',mode='w') # double acquisitions with both spacers 
    K = open(outPath+outName+'.all.fasta',mode='w')
    MC = open(outPath+outName+'.multipleAcquisitions.complete.fasta',mode='w') # multiple acquisitions with all spacers
    M = open(outPath+outName+'.multipleAcquisitions.fasta',mode='w')
    if not os.path.exists('outputs'):
        os.makedirs('outputs')
    SS = open("outputs/summaryStats.txt", mode = 'a')
    if infoFile:
        GC = open(outPath+outName+'.info.txt',mode='w')
        GC.write("Unique_Spacer_Sequence"+'\t'+"Sequence_Length"+'\t'+"GC_content"+'\n')


    os.system(str("echo '##################################################'"))
    os.system(str('echo '+"'"+inFile+' accepted for processing'+"'"))

    readName = ''
    D={}
    rawReads=0 # total reads in fasta
    spacerReads=0 # reads having a spacer
    UniqueSingleAcquisitions=0 # number of unique single acquisitions based on spacer sequence
    UniqueDoubleAcquisitions=0 # number of unique double acquisitions based on spacer sequence
    UniqueMultipleAcquisitions=0 # number of unique multiple acquisitions based on spacer sequence
    SinglefullRepeatReads=0 # reads with a double acquisition
    MultiplefullRepeatReads=0 # reads with multiple acquisitions
    spacerReadsDoubleBoth=0 # double acquisitions with both spacers
    spacerReadsDoubleOne=0 # double acquisitions with one spacer
    spacerReadsMultiComplete=0 # multiple acquisition with all spacers
    spacerReadsMultiSome=0 # multiple acquisiton with one spacer
    spacerReadsMultiNoSpacerBetweenFullDRs=0 # multiple DRs without spacers between them (probably PCR artifacts)
    minReadLength=minStagger + len(firstRepeat) + minSpacer + len(secondRepeat)
    allDistalSpacers=[]
    if 'null' not in LBC:
        NonLBCReads=0

   
    for L in F: # loop through reads in file
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            rawReads+=1
            continue
        L=L.strip()
            
        # ignore reads that are too short to detect adapted spacer
        if len(L) < minReadLength: 
            continue

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
                continue
            
        # identify and store reads with more than one acquisition (ie those that contain a full DR sequence)

        numMismatches = int(round(len(fullRepeat)*0.1)) 
        tempFullRepeat = fuzzysearch.find_near_matches(fullRepeat, L, max_l_dist = numMismatches)

        if tempFullRepeat:  # if full repeat is found

        ## DOUBLE ACQUISITIONS ##
            if len(tempFullRepeat) is 1:  

                SinglefullRepeatReads += 1
                I.write(readName+'\n'+L+'\n')
                
                # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
                tempSpacerProximal = L[:tempFullRepeat[0].start+len(secondRepeat)]  
                firstExpect = minStagger
                firstRangeToLook = maxStagger - minStagger + len(fullRepeat) - len(firstRepeat) # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                secondExpect = firstExpect + len(firstRepeat) + minSpacer
                secondRangeToLook = maxStagger - minStagger + maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                spacerProximal=editSpacer(tempSpacerProximal,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)
                
                tempSpacerDistal = L[tempFullRepeat[0].end-len(firstRepeat):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
                firstExpect = 0
                firstRangeToLook = len(fullRepeat) - len(firstRepeat) # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                secondExpect = firstExpect + len(firstRepeat) + minSpacer
                secondRangeToLook = maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                spacerDistal=editSpacer(tempSpacerDistal,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)

                # if single reads have distal & proximal spacers, label read and export these to file for looking into distal-proximal pairs

                if spacerDistal and spacerProximal:
                    spacerReads+=1
                    spacerReadsDoubleBoth+=1

                    doubleSpacer= spacerProximal+"_"+spacerDistal
                        
                    # store spacers in dict, this will force uniqueness

                    if doubleSpacer not in allDistalSpacers:
                        if doubleSpacer not in D:
                            D[doubleSpacer]=[readName+'_doubleAcquisitions',0]
                        D[doubleSpacer][1]+=1
                        allDistalSpacers.append(spacerDistal)
                        J.write(readName+'_doubleAcquisitions_both_distal'+'\n'+spacerDistal+'\n')
                        J.write(readName+'_doubleAcquisitions_both_proximal'+'\n'+spacerProximal+'\n')

                    K.write(readName+'_doubleAcquisitions_both_distal'+'\n'+spacerDistal+'\n')
                    K.write(readName+'_doubleAcquisitions_both_proximal'+'\n'+spacerProximal+'\n')

                # process reads with only distal spacer
                elif spacerDistal:
                    spacerReads+=1
                    spacerReadsDoubleOne+=1

                    if spacerDistal not in allDistalSpacers:
                        if spacerDistal not in D:
                            D[spacerDistal]=[readName+'_doubleAcquisitions_distal',0]
                        D[spacerDistal][1]+=1

                    K.write(readName+'_doubleAcquisitions_distal'+'\n'+spacerDistal+'\n')

                # process reads with only proximal spacer
                elif spacerProximal:
                    spacerReads+=1
                    spacerReadsDoubleOne+=1

                    if spacerProximal not in D:
                        D[spacerProximal]=[readName+'_doubleAcquisitions_proximal',0]
                    D[spacerProximal][1]+=1

                    K.write(readName+'_doubleAcquisitions_proximal'+'\n'+spacerProximal+'\n')
                else: 
                    continue

            ## MULTIPLE ACQUISITIONS ##

            elif len(tempFullRepeat)>1:
                MultiplefullRepeatReads+=1
                M.write(readName+'\n'+L+'\n')

                tempSpacerProximal = L[:tempFullRepeat[0].start+len(secondRepeat)]
                firstExpect = minStagger
                firstRangeToLook = maxStagger - minStagger + len(fullRepeat) - len(firstRepeat) # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                secondExpect = firstExpect + len(firstRepeat) + minSpacer
                secondRangeToLook = maxStagger - minStagger + maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                spacerProximal=editSpacer(tempSpacerProximal,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)

                spacersMedial=''

                for index, DR in enumerate(tempFullRepeat):
                    length = tempFullRepeat[index].start - tempFullRepeat[index-1].end
                    if length > minSpacer and length < maxSpacer:
                        tempspacersMedial = L[tempFullRepeat[index-1].end:tempFullRepeat[index].start]
                        spacersMedial=spacersMedial+"_"+tempspacersMedial
                if spacersMedial:
                    spacersMedial=spacersMedial[1:]
                        

                tempSpacerDistal = L[tempFullRepeat[-1].end-len(firstRepeat):]
                firstExpect = 0
                firstRangeToLook = len(fullRepeat) - len(firstRepeat) # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                secondExpect = firstExpect + len(firstRepeat) + minSpacer
                secondRangeToLook = maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                spacerDistal=editSpacer(tempSpacerDistal,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)

                if spacerDistal and spacerProximal and spacersMedial:

                    multipleSpacer=spacerProximal+"_"+spacersMedial+"_"+spacerDistal
                    spacerReads+=1
                    spacerReadsMultiComplete+=1

                    if multipleSpacer not in allDistalSpacers:

                        if multipleSpacer not in D:
                            D[multipleSpacer]=[readName+'_multipleAcquisitions',0]
                        D[multipleSpacer][1]+=1

                        for i in range(1,multipleSpacer.count('_')+1):
                            if multipleSpacer.split('_', i)[i] not in allDistalSpacers:
                                allDistalSpacers.append(multipleSpacer.split('_', i)[i])

                        MC.write(readName+'_multipleAcquisitions_complete_distal'+'\n'+spacerDistal+'\n')
                        MC.write(readName+'_multipleAcquisitions_complete_proximal'+'\n'+spacerProximal+'\n')
                        for index,spacer in enumerate(spacersMedial.split('_')):
                            MC.write(readName+'_multipleAcquisitions_complete_medial_'+str(index)+'\n'+spacer+'\n')

                    K.write(readName+'_multipleAcquisitions_complete_distal'+'\n'+spacerDistal+'\n')
                    K.write(readName+'_multipleAcquisitions_complete_proximal'+'\n'+spacerProximal+'\n')

                    for index,spacer in enumerate(spacersMedial.split('_')):
                        K.write(readName+'_multipleAcquisitions_complete_medial_'+str(index)+'\n'+spacer+'\n')

                elif spacerDistal and spacersMedial:

                    partMultipleSpacer = spacersMedial + "_" + spacerDistal
                    spacerReads+=1
                    spacerReadsMultiSome+=1
                    
                    if partMultipleSpacer not in allDistalSpacers:
                        if partMultipleSpacer not in D:
                            D[partMultipleSpacer]=[readName+'_multipleAcquisitions_noProximal',0]
                        D[partMultipleSpacer][1]+=1

                        for i in range(1,partMultipleSpacer.count('_')+1):
                            if partMultipleSpacer.split('_', i)[i] not in allDistalSpacers:
                                allDistalSpacers.append(partMultipleSpacer.split('_', i)[i])


                    K.write(readName+'_multipleAcquisitions_noProximal_distal'+'\n'+spacerDistal+'\n')
                    for index,spacer in enumerate(spacersMedial.split('_')):
                        K.write(readName+'_multipleAcquisitions_noProximal_medial_'+str(index)+'\n'+spacer+'\n')

                elif spacerProximal and spacersMedial:
                    spacerReads+=1
                    spacerReadsMultiSome+=1
                    partMultipleSpacer =spacerProximal+"_"+spacersMedial

                    if partMultipleSpacer not in allDistalSpacers:
                        if partMultipleSpacer not in D:
                            D[partMultipleSpacer]=[readName+'_multipleAcquisitions_noDistal',0]
                        D[partMultipleSpacer][1]+=1

                        for i in range(1,partMultipleSpacer.count('_')+1):
                            if partMultipleSpacer.split('_', i)[i] not in allDistalSpacers:
                                allDistalSpacers.append(partMultipleSpacer.split('_', i)[i])

                    K.write(readName+'_multipleAcquisitions_noDistal_proximal'+'\n'+spacerProximal+'\n')
                    for index,spacer in enumerate(spacersMedial.split('_')):
                        K.write(readName+'_multipleAcquisitions_noDistal_medial_'+str(index)+'\n'+spacer+'\n')

                elif spacerDistal and spacerProximal and not spacersMedial:
                    spacerReads+=1
                    spacerReadsMultiSome+=1
                    spacerReadsMultiNoSpacerBetweenFullDRs+=1
                    partMultipleSpacer=spacerProximal+"_"+spacerDistal

                    if partMultipleSpacer not in allDistalSpacers:
                        if partMultipleSpacer not in D:
                            D[partMultipleSpacer]=[readName+'_multipleAcquisitions_noMedial',0]
                        D[partMultipleSpacer][1]+=1
                        if spacerDistal not in allDistalSpacers:
                            allDistalSpacers.append(spacerDistal)

                    K.write(readName+'_multipleAcquisitions_noMedial_distal'+'\n'+spacerDistal+'\n')
                    K.write(readName+'_multipleAcquisitions_noMedial_proximal'+'\n'+spacerProximal+'\n')

                elif spacerDistal:
                    spacerReads+=1
                    spacerReadsMultiSome+=1
                    spacerReadsMultiNoSpacerBetweenFullDRs+=1
                    if spacerDistal not in allDistalSpacers:
                        if spacerDistal not in D:
                            # store spacers in dict, this will force uniqueness
                            D[spacerDistal]=[readName+'_multipleAcquisitions_onlyDistal',0]
                        D[spacerDistal][1]+=1

                    K.write(readName+'_multipleAcquisitions_onlyDistal'+'\n'+spacerDistal+'\n')

                elif spacerProximal:
                    spacerReads+=1
                    spacerReadsMultiSome+=1
                    spacerReadsMultiNoSpacerBetweenFullDRs+=1

                    if spacerProximal not in D:
                        # store spacers in dict, this will force uniqueness
                        D[spacerProximal]=[readName+'_multipleAcquisitions_onlyProximal',0]
                    D[spacerProximal][1]+=1

                    K.write(readName+'_multipleAcquisitions_onlyProximal'+'\n'+spacerProximal+'\n')

                elif spacersMedial:
                    spacerReads+=1
                    spacerReadsMultiSome+=1

                    partMultipleSpacer = spacersMedial

                    if partMultipleSpacer not in allDistalSpacers:
                        if partMultipleSpacer not in D:
                            D[partMultipleSpacer]=[readName+'_multipleAcquisitions_onlyMedial',0]
                        D[partMultipleSpacer][1]+=1

                        for i in range(1,partMultipleSpacer.count('_')+1):
                            if partMultipleSpacer.split('_', i)[i] not in allDistalSpacers:
                                allDistalSpacers.append(partMultipleSpacer.split('_', i)[i])

                    for index,spacer in enumerate(spacersMedial):
                        K.write(readName+'_multipleAcquisitions_only_medial_'+str(index)+'\n'+spacer+'\n')
            else:
                continue

        # run standard code if multiple acquisitions (ie full DR sequence) not detected
        else:
            firstExpect = minStagger
            firstRangeToLook = maxStagger - minStagger + len(fullRepeat) - len(firstRepeat)  # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
            secondExpect = firstExpect + len(firstRepeat) + minSpacer
            secondRangeToLook = maxStagger - minStagger + maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
            spacer=editSpacer(L,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)
            if not spacer: 
                continue
            spacerReads+=1
            if spacer not in allDistalSpacers:
                if spacer not in D:
                    D[spacer]=[readName,0]
                D[spacer][1]+=1

            K.write(readName+'\n'+spacer+'\n')

    # iterate through spacers and print to file
    for spacer in allDistalSpacers:
        if spacer in D:
            del(D[spacer])

    seqlist=sorted(D.keys())    
    for S in seqlist:
        if 'double' in D[S][0]:
            UniqueDoubleAcquisitions+=1
        elif 'multiple' in D[S][0]:
            UniqueMultipleAcquisitions+=1
        else:
            UniqueSingleAcquisitions+=1
        if '_' not in S and len(S)>0:
            G.write(D[S][0]+'_rep_'+str(D[S][1])+'\n'+S+'\n')
            if infoFile:
                count = Counter(S)
                gc = (count['G']+count['C'])*100/len(S)
                GC.write(str(S)+'\t'+str(len(S))+'\t'+str(gc)+'\n')

        else:
            for index, spcr in enumerate(S.split('_')):
                if len(spcr)>0:
                    G.write(D[S][0]+'_rep_'+str(D[S][1])+'_spacerPosition_'+str(index)+'\n'+spcr+'\n')
                    if infoFile:
                        count = Counter(spcr)
                        gc = (count['G']+count['C'])*100/len(spcr)
                        GC.write(str(spcr)+'\t'+str(len(spcr))+'\t'+str(gc)+'\n')


    if not D:
        os.system(str("echo 'No spacers to map'"))

    NonLBCReadPercentage = 'null'
    if 'null' not in LBC:
        NonLBCReadPercentage = NonLBCReads*100/rawReads
    os.system(str('echo '+"'"+'*'+'\t'+' sampleName'+'\t'+' rawReads'+'\t'+'nonLibraryBarcodeRead%'+'\t'+'identifiedSpacers'+'\t'+'uniqueSpacers'+'\t'+'doubleAcquisitions'+'\t'+'doubleAcquisitions.paired'+'\t'+'multipleAcquisitions'+'\t'+'multipleAcquisitions.complete'+'\t'+'uniqueSingleAcquisitions'+'\t'+'uniqueDoubleAcquisitions'+'\t'+'uniqueMultipleAcquisitions'+"'"))
    os.system(str('echo '+"'"+'@'+'\t'+str(outPath+outName)+'\t'+str(rawReads)+'\t'+str(NonLBCReadPercentage)+'\t'+str(spacerReads)+'\t'+str(len(D.keys()))+'\t'+str(SinglefullRepeatReads)+'\t'+str(spacerReadsDoubleBoth)+'\t'+str(MultiplefullRepeatReads)+'\t'+str(spacerReadsMultiComplete)+'\t'+str(UniqueSingleAcquisitions)+'\t'+str(UniqueDoubleAcquisitions)+'\t'+str(UniqueMultipleAcquisitions)+"'"))
    SS.write(str(outName)+'\t'+str(rawReads)+'\t'+str(NonLBCReadPercentage)+'\t'+str(spacerReads)+'\t'+str(len(D.keys()))+'\t'+str(SinglefullRepeatReads)+'\t'+str(spacerReadsDoubleBoth)+'\t'+str(MultiplefullRepeatReads)+'\t'+str(spacerReadsMultiComplete)+'\t'+str(UniqueSingleAcquisitions)+'\t'+str(UniqueDoubleAcquisitions)+'\t'+str(UniqueMultipleAcquisitions)+'\n')
   

      
    # print useful summary stats to stdout

    F.close()
    G.close()
    I.close()
    J.close()
    K.close()
    MC.close()
    M.close()
    SS.close()
    if infoFile:
        GC.close()
    