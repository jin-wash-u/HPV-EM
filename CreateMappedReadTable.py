#!/usr/bin/env python
#Create table of mapped read matches and mismatches for use as input to HPV type EM algorithm
from __future__ import print_function
import sys
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE

#Calculate score to identify low-complexity reads using DUST algorithm
#(S>2 should be filtered)
def dust(read):
    tripletDict = {}
    for i in range(len(read)-2):
        c = read[i:i+3]
        if c in tripletDict:
            tripletDict[c] += 1
        else:
            tripletDict[c] = 1
    S = 0
    l = len(read)-2
    for trip in tripletDict:
        c = float(tripletDict[trip])
        S += c*(c-1)/2/(l-1)
    return S


def mapReads(hpvBams, defaultHpvRef=True, hpvRefPath='', filterLowComplex=True, outputName='hpvType'):
    mapped_reads = set()
    hpvTypeDict = {}
    hpvRefIdDict = {}
    hpvRefSeqDict = {}
    hpvRefCovDict = {}
    installDir = os.path.dirname(os.path.abspath(__file__))

    #Make dict to translate ref seq names (SAM field 2) into HPV type names
    if defaultHpvRef:
        hpvRefPath = installDir+'/reference/combined_pave_hpv.fa'
        with open(installDir+'/reference/hpv_type_dict.tsv','r') as fHpvNames:
            for line in fHpvNames:
                line = line.strip().split('\t')
                hpvRefIdDict[line[0]] = line[1]

    #Read in HPV reference file
    with open(hpvRefPath,'r') as fHpvRef:
        hpvRef = ''
        for line in fHpvRef:
            if not line:
                break
            
            if line[0]=='>':
                if hpvRef:
                    hpvRefSeqDict[refId] = hpvRef
                    hpvRef = ''
                refId = line.strip().split()[0][1:]
            else:
                hpvRef+=line.strip()
        hpvRefSeqDict[refId] = hpvRef

    #For all HPV*.sam files in directory
    for bam in hpvBams:
        mate = bam.split('.')[2]

        ##Read the file
        cmd = ['samtools','view', bam]
        pipe = Popen(cmd, stdout=PIPE)
        ##loop over lines
        for line in pipe.stdout:
            #for each line, check if field 5 == 90M (or [readLen]M)
            #if so, get name from field 0, ref id from field 2, seq from field 9, and position from field 3
            line = line.strip().split('\t')
            [readName,refIdList,readPosList,readCIGARList,readSeq] = [line[0],[line[2]],[int(line[3])],[line[5]],line[9]]
            readLen = len(readSeq)
            readSeq = readSeq.upper()
            readName+='/'+mate

            if dust(readSeq)<=2 or not filterLowComplex:
                #Get any alternative alignments
                if line[-1][:2]=='XA':
                    altAligns = line[-1][5:].split(';')[:-1]
                    for align in altAligns:
                        align = align.split(',')
                        refIdList.append(align[0])
                        readPosList.append(int(align[1][1:]))
                        readCIGARList.append(align[2])
                        #Last field [3] is edit distance (as in NM tag)
                for i in range(len(refIdList)):
                    refId = refIdList[i]
                    readPos = readPosList[i]
                    readCIGAR = readCIGARList[i]
                    if readCIGAR == str(readLen)+'M':
                        #Fetch hpvRef from dict if already loaded, otherwise read it in
                        if refId in hpvRefSeqDict:
                            hpvRef = hpvRefSeqDict[refId]
                        else:
                            raise KeyError('No reference genome available for genotype '+refId+', please check genome names match BAM file')

                        #Add name to set of mapped_reads names
                        mapped_reads.add(readName)
                        if refId not in hpvRefCovDict:
                            hpvRefCovDict[refId] = [0]*len(hpvRefSeqDict[refId])
                        hpvRefCovDict[refId][(readPos-1):(readPos-1+readLen)] = [c+1 for c in hpvRefCovDict[refId][(readPos-1):(readPos-1+readLen)]]
                        refSeq = hpvRef[(readPos-1):(readPos-1+readLen)].upper()
                        Le = 0
                        for j in range(readLen):
                            if readSeq[j]!=refSeq[j]:
                                Le+=1
                        Lm = readLen-Le

                        #If not in there, add dict for this HPV type to dict of dicts
                        if refId not in hpvTypeDict:
                            hpvTypeDict[refId] = {}
                        #Add read to dict of reads for this HPV type, value is [Lm,Le]
                        hpvTypeDict[refId][readName]=[Lm,Le]
        pipe.stdout.close()

    totalReads = len(mapped_reads)
    #Write out table for EM algo
    outTable = [totalReads]
    for hpv in hpvTypeDict:
        if defaultHpvRef:
            hpvType = hpvRefIdDict[hpv].replace('REF','').replace('.fa','')
        else:
            hpvType = hpv
        outLine = hpvType
        for read in mapped_reads:
            if read in hpvTypeDict[hpv]:
                Lm = hpvTypeDict[hpv][read][0]
                Le = hpvTypeDict[hpv][read][1]
                outLine += '\t1\t'+str(Lm)+'\t'+str(Le)
            else:
                outLine += '\t0\t-1\t-1'
        outTable.append(outLine)

    #Plot coverage maps
    for hpv in hpvRefCovDict:
        fig = plt.figure()
        cov = fig.add_subplot(111)
        cov.plot(list(range(len(hpvRefCovDict[hpv]))), hpvRefCovDict[hpv], 'k', lw = 0.8)
        cov.set_xlabel('Nucleotide', fontsize = 14, color = 'black')
        cov.set_ylabel('Read coverage', fontsize = 14, color = 'black')
        if defaultHpvRef:
            hpvName = hpvRefIdDict[hpv].replace('REF','').replace('.fa','')
        else:
            hpvName = hpv.replace(' ','')
        plt.title(hpvName)
        fig.tight_layout()
        fig.savefig(outputName+'.'+hpvName+'.cov.pdf',bbox_inches='tight')
        plt.close(fig)
    return outTable


def main(argv):
    outTable = mapReads(argv[1:3])
    with open('hpvType.mappedReads.tsv','w') as outFile:
        for line in outTable:
            outFile.write(str(line)+'\n')

        
if __name__=="__main__":
    main(sys.argv)
