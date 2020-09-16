#!/usr/bin/env python
#Create table of mapped read matches and mismatches for use as input to HPV type EM algorithm
import sys
import os
import re
import time
import argparse as argp
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, lines as lines
from subprocess import Popen, PIPE

class alignInfo:
    def __init__(self, Lm, Le, pos, cigar):
        self.Lm = Lm
        self.Le = Le
        self.pos = pos
        self.cigar = cigar
        
class readAligns:
    def __init__(self, refId, passDust, mate, Lm, Le, pos, cigar):
        mate = int(mate) - 1
        Lm = int(Lm)
        Le = int(Le)
        pos = int(pos)
        self.isAmbig = False
        self.passDust = [False, False]
        self.passDust[mate] = passDust
        self.dictRefId_AlignInfo = {}
        self.dictRefId_AlignInfo[refId] = [0,0]
        self.dictRefId_AlignInfo[refId][mate] = alignInfo(Lm, Le, pos, cigar)

    def addAlign(self, refId, passDust, mate, Lm, Le, pos, cigar):
        mate = int(mate) - 1
        Lm = int(Lm)
        Le = int(Le)
        pos = int(pos)
        if refId in self.dictRefId_AlignInfo:
            if self.dictRefId_AlignInfo[refId][mate]:
                if Lm > self.dictRefId_AlignInfo[refId][mate].Lm:
                    self.dictRefId_AlignInfo[refId][mate] = alignInfo(Lm, Le, pos, cigar)
                    self.passDust[mate] = self.passDust[mate] or passDust
            else:
                self.dictRefId_AlignInfo[refId][mate] = alignInfo(Lm, Le, pos, cigar)
                self.passDust[mate] = self.passDust[mate] or passDust
        else:
            self.isAmbig = True
            self.passDust[mate] = self.passDust[mate] or passDust
            self.dictRefId_AlignInfo[refId] = [0,0]
            self.dictRefId_AlignInfo[refId][mate] = alignInfo(Lm, Le, pos, cigar)


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


def mapReads(hpvBams, defaultHpvRef=True, hpvRefPath='', filterLowComplex=True, outputName='hpvType', covMapYmax=0):
    mapped_reads = set()
    dictReadName_ReadAligns = {}
    hpvRefIdMappedSet = set()
    hpvRefIdMappedNumDict = {}
    hpvRefIdGeneDict = {}
    hpvRefIdSeqDict = {}
    hpvRefIdCovDict = {}

    installDir = os.path.dirname(os.path.abspath(__file__))

    # Make dict to translate ref seq names (SAM field 2) into HPV type names
    if defaultHpvRef:
        hpvRefPath = installDir+'/reference/combined_pave_hpv.fa'
        
    if defaultHpvRef:
        with open(installDir+'/reference/hpv_gene_annot.tsv','r') as fHpvGenes:
            for line in fHpvGenes:
                line = line.strip().split('\t')
                if line[0] in hpvRefIdGeneDict:
                    hpvRefIdGeneDict[line[0]].append(line[1:])
                else:
                    hpvRefIdGeneDict[line[0]] = [line[1:]]

    annotColorDict = {'E1':'g','E2':'gray','E3':'y','E4':'r','E5':'orange',
                      'E6':'b','E7':'m','E8':'c','L1':'indigo','L2':'brown'}
    annotColors=['maroon','navy','pink','g','gray','k','y','r','orange','b','m','c','indigo']

    # Read in HPV reference file
    with open(hpvRefPath,'r') as fHpvRef:
        hpvRef = ''
        refId = ''
        for line in fHpvRef:
            if not line:
                break
            
            if line[0]=='>':
                if hpvRef:
                    hpvRefIdSeqDict[refId] = hpvRef
                    hpvRef = ''
                refId = line.strip().split()[0][1:]
            else:
                hpvRef+=line.strip()
        hpvRefIdSeqDict[refId] = hpvRef

    # For all HPV*.bam files in directory
    for bam in hpvBams:
        mate = bam.split('.')[-4]
        
        # Read the file
        cmdArgs = ['samtools','view', bam]
        if sys.version[0] == '2':
            pipe = Popen(cmdArgs, stdout=PIPE)
        else:
            pipe = Popen(cmdArgs, stdout=PIPE, encoding='utf8')
        # loop over lines
        for line in pipe.stdout:
            # Get read name from field 0, SAM flags from f1, ref id from f2,
            # position from f3, seq from field 9, and tags from field 11
            line = line.strip().split('\t')
            [readName,readFlags,readRefId,readPos,readCIGAR,readSeq,readTags] = \
                     [line[0],line[1],line[2],int(line[3]),line[5],line[9],line[11:]]
            readLen = len(readSeq)
            readSeq = readSeq.upper()
            try:
                editDist = [tag for tag in readTags if tag.startswith('NM')][0].split(':')[-1]
            except:
                print('Error parsing tags:')
                print(readName)
                print(readTags)
                sys.exit(1)
            
            # Add this read to the dictionary
            cigarList = list(filter(None, re.split('(\D+)',readCIGAR)))
            alignedSeq = ''
            pos=0
            clipLen=0
            for cigar in zip(cigarList[0::2], cigarList[1::2]):
                clen = int(cigar[0])
                if cigar[1] in 'M=XIP':
                    alignedSeq += readSeq[pos:pos+clen]
                    pos += clen
                elif cigar[1] in 'SH':
                    pos += clen

            # Get proper length of matching using corrected readlength
            Le = int(editDist)
            Lm = len(alignedSeq) - Le

            passDust = dust(alignedSeq)<=2
            # Disallow clipping on both ends
            if cigarList[1] in 'HS' and cigarList[-1] in 'HS':
                passDust = False
            hpvRefIdMappedSet.add(readRefId)
            if readName in dictReadName_ReadAligns:
                dictReadName_ReadAligns[readName].addAlign(readRefId, passDust, mate, Lm, Le, readPos, readCIGAR)
            else:
                dictReadName_ReadAligns[readName] = readAligns(readRefId, passDust, mate, Lm, Le, readPos, readCIGAR)
        while pipe.poll() is None:
            # Process not yet terminated, wait
            time.sleep(0.5)
        if pipe.returncode > 0:
            raise RuntimeError('Error parsing viral-aligned BAM files; aborting.')
        
    # Get numbers of transcripts mapped to each type
    for refId in hpvRefIdMappedSet:
        hpvRefIdMappedNumDict[refId] = 0
    for readName in dictReadName_ReadAligns.keys():
        ra = dictReadName_ReadAligns[readName]
        for refId in hpvRefIdMappedSet:
            if refId in ra.dictRefId_AlignInfo:
                hpvRefIdMappedNumDict[refId] += 1
                
    # Check if all reads aligned to an HPV type have equal or better alignment to a type with more reads
    refIdList = sorted(hpvRefIdMappedSet)
    for refId in refIdList:
        isRedundant = True
        for readName in dictReadName_ReadAligns.keys():
            ra = dictReadName_ReadAligns[readName]
            if refId in ra.dictRefId_AlignInfo:
                if not ra.isAmbig:
                    isRedundant = False
                    break
                
        if isRedundant:
            for readName in dictReadName_ReadAligns.keys():
                ra = dictReadName_ReadAligns[readName]
                if refId in ra.dictRefId_AlignInfo:
                    rai = ra.dictRefId_AlignInfo
                    thisLm = 0
                    if rai[refId][0]:
                        thisLm += rai[refId][0].Lm
                    if rai[refId][1]:
                        thisLm += rai[refId][1].Lm
                    anyRedundant = False
                    for altRefId in rai:
                        altLm = 0
                        if rai[altRefId][0]:
                            altLm += rai[altRefId][0].Lm
                        if rai[altRefId][1]:
                            altLm += rai[altRefId][1].Lm
                        if altLm >= thisLm and hpvRefIdMappedNumDict[altRefId] > hpvRefIdMappedNumDict[refId]:
                            anyRedundant = True
                    if anyRedundant:
                        del rai[refId]
                        hpvRefIdMappedNumDict[refId] -= 1

    #Now process all reads/read pairs in dict to prepare output table and coverage maps
    # Each column of the outTable is a distinct, mapped read
    # First line of outTable is total mapped read #, followed by unique(U)/ambiguous(A) status of each read/pair
    # There follows 1 line for each HPV reference with at least one read mapped to it.  The first column is the HPV name,
    # and each read column has the following format: [0/1 (whether maps to this reference), Lm (-1 if unmapped), ...
    #  ... Le (-1 if unmapped), (comma-separated gene list)]
    # First line output
    mappedCount = 0
    outLine = ''
    nameLine = ''
    for readName in list(dictReadName_ReadAligns.keys()):
        ra = dictReadName_ReadAligns[readName]
        if filterLowComplex and not any(ra.passDust):
            del dictReadName_ReadAligns[readName]
        else:
            mappedCount += 1
            if ra.isAmbig:
                outLine += '\tA'
            else:
                outLine += '\tU'
            nameLine += '\t'+readName
    outLine = str(mappedCount)+outLine
    outTable = [nameLine]
    outTable.append(outLine)

    # Rest of table
    if mappedCount:
        for refId in hpvRefIdMappedSet:
            if defaultHpvRef:
                if refId.split('.')[0][-3:]=='REF':
                    hpvName = refId.split('.')[0][:-3]
                elif refId.split('.')[0][-2:]=='nr':
                    hpvName = refId.split('.')[0][:-2]
                else:
                    hpvName = refId.replace(' ','')
            else:
                hpvName = refId.replace(' ','')    
            outLine = hpvName
            for readName in dictReadName_ReadAligns:
                ra = dictReadName_ReadAligns[readName]
                if refId in ra.dictRefId_AlignInfo:
                    geneSet = set()
                    Lm = 0
                    Le = 0

                    # Update read coverage depths for this HPV type
                    if refId not in hpvRefIdCovDict:
                        hpvRefIdCovDict[refId] = [0]*len(hpvRefIdSeqDict[refId])
                    for mInd in range(2):
                        mate = ra.dictRefId_AlignInfo[refId][mInd]
                        if mate:
                            cigarList = list(filter(None, re.split('(\D+)',mate.cigar)))
                            pos = mate.pos
                            for cigar in zip(cigarList[0::2], cigarList[1::2]):
                                if cigar[1] in 'M=X':
                                    for i in range(int(cigar[0])):
                                        # Add to coverage count
                                        try:
                                            hpvRefIdCovDict[refId][pos-1] += 1
                                        except:
                                            print('readName: {}; refID: {}; startPos: {}; CIGAR: {}; pos: {}'.format(readName,refId,mate.pos,mate.cigar,pos))
                                            print('Len(hpvRefIdCovDict[refId]): {}'.format(len(hpvRefIdCovDict[refId])))
                                            raise

                                        # Mark any genes this read covers
                                        for gene in hpvRefIdGeneDict[refId]:
                                            gName = gene[0]
                                            gStart = int(gene[1])
                                            gEnd = int(gene[2])
                                            if gStart <= pos and pos <= gEnd:
                                                geneSet.add(gName)

                                        pos = pos+1
                                elif cigar[1] in 'DN':
                                    for i in range(int(cigar[0])):
                                        pos = pos+1
                                #else :'IPSH'

                            Lm += mate.Lm
                            Le += mate.Le
                    genes = ','.join(sorted(geneSet))
                    outLine += '\t'+'\t'.join(['1',str(Lm),str(Le),genes])
                else:
                    outLine += '\t0\t-1\t-1\t'
            outTable.append(outLine)

    # Plot coverage maps
    for refId in hpvRefIdCovDict:
        fig = plt.figure(figsize=(9,4))
        r=fig.canvas.get_renderer()
        cov = fig.add_subplot(111)
        cov.plot(list(range(len(hpvRefIdCovDict[refId]))), hpvRefIdCovDict[refId], 'k', lw = 0.8)
        cov.set_ylabel('Read coverage', fontsize = 14, color = 'black')
        if defaultHpvRef:
            if refId.split('.')[0][-3:]=='REF':
                hpvName = refId.split('.')[0][:-3]
            elif refId.split('.')[0][-2:]=='nr':
                hpvName = refId.split('.')[0][:-2]
            else:
                hpvName = refId.replace(' ','')
        else:
            hpvName = refId.replace(' ','')
        plt.title(hpvName)

        if covMapYmax:
            cov.set_ylim(top=covMapYmax)

        # Plot gene annotations
        glines = []
        glabels = []
        y1end=0
        y2end=0
        annotScale = 1.3
        ypos1 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/12*annotScale
        ypos2 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/7.9*annotScale
        ypos3 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/5.8*annotScale
        yposlab1 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/8.5*annotScale
        yposlab2 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/6.2*annotScale
        yposlab3 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/4.8*annotScale
        if refId in hpvRefIdGeneDict:
            ic = 0
            gNameLast = ''
            for gene in hpvRefIdGeneDict[refId]:
                gName = gene[0]
                gStart = int(gene[1])
                gEnd = int(gene[2])
                
                tname1 = gName[:2].upper()
                tname2 = gName[-2:].upper()
                if (tname1 in annotColorDict and
                    (len(gName)<3 or gName[2] not in '^*')):
                    gc = annotColorDict[tname1]
                elif (tname2 in annotColorDict and
                      (len(gName)<3 or gName[-3] not in '^*')):
                    gc = annotColorDict[tname2]
                else:
                    if gName != gNameLast:
                        ic += 1
                    gc = annotColors[ic]
                    if ic>13:
                        ic = 0
                if gStart >= y1end:
                    ypos = ypos1
                    yposlab = yposlab1
                elif gStart >= y2end:
                    ypos = ypos2
                    yposlab = yposlab2
                else:
                    ypos = ypos3
                    yposlab = yposlab3
                gline = cov.add_line(lines.Line2D([gStart,gEnd],[ypos,ypos],color=gc,clip_on=False, linewidth=2))
                glines.append(gline)
                glabel = cov.text(gStart, yposlab, gName)
                glabels.append(glabel)

                if ypos == ypos1:
                    y1end = max(gEnd,
                                cov.transData.inverted().transform(glabel.get_window_extent(renderer=r))[1][0])
                elif ypos == ypos2:
                    y2end = max(gEnd,
                                cov.transData.inverted().transform(glabel.get_window_extent(renderer=r))[1][0])
                gNameLast = gName

        fig.savefig(outputName+'.'+hpvName+'.cov.pdf',bbox_inches='tight',bbox_extra_artists=glines+glabels)
        plt.close(fig)
    return outTable


def main(argv):
    mapParse = argp.ArgumentParser()
    mapParse.add_argument('bam1')
    mapParse.add_argument('bam2', nargs='?', help='(optional)', default='not supplied')
    mapParse.add_argument('-r','--reference', default=0)
    mapParse.add_argument('-o','--outname', type=str, default='./hpvType')
    mapParse.add_argument('-d','--disabledust', action='store_true')
    mapParse.add_argument('-y','--ylimit', type=int, help='fix a maximum y-value for all coverage map axes', default=0) 
    args = mapParse.parse_args()

    hpvBams = [args.bam1]
    if args.bam2 != "not supplied":
        hpvBams += [args.bam2]
        
    if(args.reference == 0):
        defaultHpvRef = True
    else:
        defaultHpvRef = False
        
    outTable = mapReads(hpvBams, defaultHpvRef=defaultHpvRef, hpvRefPath=args.reference, filterLowComplex=not(args.disabledust), outputName=args.outname, covMapYmax=args.ylimit)
    with open(args.outname+'.mappedReads.tsv','w') as outFile:
        for line in outTable:
            outFile.write(str(line)+'\n')

        
if __name__=="__main__":
    main(sys.argv)
