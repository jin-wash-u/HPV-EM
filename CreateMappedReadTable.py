#!/usr/bin/env python
#Create table of mapped read matches and mismatches for use as input to HPV type EM algorithm
#from __future__ import print_function
import sys
import os
import re
import time
import argparse as argp
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, lines as lines #import matplotlib.pyplot as plt
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


def mapReads(hpvBams, defaultHpvRef=True, hpvRefPath='', filterLowComplex=True, outputName='hpvType', covMapYmax=0):
    mapped_reads = set()
    hpvTypeDict = {}
    hpvRefIdDict = {}
    hpvGeneDict = {}
    hpvRefSeqDict = {}
    hpvRefCovDict = {}
    hpvRefReadCountsDict = {}
    installDir = os.path.dirname(os.path.abspath(__file__))

    #Make dict to translate ref seq names (SAM field 2) into HPV type names
    if defaultHpvRef:
        hpvRefPath = installDir+'/reference/combined_pave_hpv.fa'
        with open(installDir+'/reference/hpv_type_dict.tsv','r') as fHpvNames:
            for line in fHpvNames:
                line = line.strip().split('\t')
                hpvRefIdDict[line[0]] = line[1]

    if defaultHpvRef:
        with open(installDir+'/reference/hpv_gene_annot.tsv','r') as fHpvGenes:
            for line in fHpvGenes:
                line = line.strip().split('\t')
                if line[0] in hpvGeneDict:
                    hpvGeneDict[line[0]].append(line[1:])
                else:
                    hpvGeneDict[line[0]] = [line[1:]]

    annotColorDict = {'E1':'g','E2':'gray','E3':'y','E4':'r','E5':'orange',
                      'E6':'b','E7':'m','E8':'c','L1':'indigo','L2':'brown'}
    annotColors=['maroon','navy','pink','g','gray','k','y','r','orange','b','m','c','indigo']

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

    #For all HPV*.bam files in directory
    for bam in hpvBams:
        mate = bam.split('.')[-2]
        
        ##Read the file
        cmdArgs = ['samtools','view', bam]
        pipe = Popen(cmdArgs, stdout=PIPE)
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

                        #Update read coverage depths for this HPV type
                        if refId not in hpvRefCovDict:
                            hpvRefCovDict[refId] = [0]*len(hpvRefSeqDict[refId])
                        hpvRefCovDict[refId][(readPos-1):(readPos-1+readLen)] = [c+1 for c in hpvRefCovDict[refId][(readPos-1):(readPos-1+readLen)]]

                        #Update gene counts for this HPV type
                        if defaultHpvRef:
                            if hpvRefIdDict[refId][-6:]=='REF.fa':
                                hpvName = hpvRefIdDict[refId][:-6]
                            elif hpvRefIdDict[refId][-5:]=='nr.fa':
                                hpvName = hpvRefIdDict[refId][:-5]
                            else:
                                hpvName = refId.replace(' ','')
                        else:
                            hpvName = refId.replace(' ','')

                        if hpvName not in hpvRefReadCountsDict:
                            hpvRefReadCountsDict[hpvName] = {}
                            if hpvName in hpvGeneDict:
                                for gene in hpvGeneDict[hpvName]:
                                    gName = gene[0]
                                    hpvRefReadCountsDict[hpvName][gName] = 0
                                    
                        for gene in hpvGeneDict[hpvName]:
                            gName = gene[0]
                            gStart = int(gene[1])
                            gEnd = int(gene[2])
                            rStart = readPos
                            rEnd = readPos+readLen
                            if gStart <= rEnd and rStart <= gEnd:
                                hpvRefReadCountsDict[hpvName][gName] += 1
                                
                        #Compare read to reference sequence, get number of matching (Lm) and mismatched (Le) bases
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
        while pipe.poll() is None:
            #Process not yet terminated, wait
            time.sleep(0.5)
        if pipe.returncode > 0:
            raise RuntimeError('Error parsing viral-aligned BAM files; aborting.')
        #pipe.stdout.close()

    totalReads = len(mapped_reads)
    #Write out table for EM algo
    outTable = [totalReads]
    for hpv in hpvTypeDict:
        if defaultHpvRef:
            if hpvRefIdDict[hpv][-6:]=='REF.fa':
                hpvName = hpvRefIdDict[hpv][:-6]
            elif hpvRefIdDict[hpv][-5:]=='nr.fa':
                hpvName = hpvRefIdDict[hpv][:-5]
            else:
                hpvName = hpv
        else:
            hpvName = hpv
        outLine = hpvName
        for read in mapped_reads:
            if read in hpvTypeDict[hpv]:
                Lm = hpvTypeDict[hpv][read][0]
                Le = hpvTypeDict[hpv][read][1]
                outLine += '\t1\t'+str(Lm)+'\t'+str(Le)
            else:
                outLine += '\t0\t-1\t-1'
        outTable.append(outLine)

    #Output read counts tables
    with open(outputName+'.readCounts.tsv','w') as fCounts:
        for hpvName in sorted(hpvRefReadCountsDict):
            fCounts.write(hpvName+'\n')
            outline1=''
            outline2=''
            for gene in sorted(hpvRefReadCountsDict[hpvName]):
                outline1 += gene+'\t'
                outline2 += str(hpvRefReadCountsDict[hpvName][gene])+'\t'
            outline1 = outline1[:-1]+'\n'
            outline2 = outline2[:-1]+'\n'
            fCounts.write(outline1+outline2+'\n')

    #Plot coverage maps
    for hpv in hpvRefCovDict:
        fig = plt.figure()
        r=fig.canvas.get_renderer()
        cov = fig.add_subplot(111)
        cov.plot(list(range(len(hpvRefCovDict[hpv]))), hpvRefCovDict[hpv], 'k', lw = 0.8)
        #cov.set_xlabel('Nucleotide', fontsize = 14, color = 'black')
        cov.set_ylabel('Read coverage', fontsize = 14, color = 'black')
        if defaultHpvRef:
            if hpvRefIdDict[hpv][-6:]=='REF.fa':
                hpvName = hpvRefIdDict[hpv][:-6]
            elif hpvRefIdDict[hpv][-5:]=='nr.fa':
                hpvName = hpvRefIdDict[hpv][:-5]
            else:
                hpvName = hpv.replace(' ','')
        else:
            hpvName = hpv.replace(' ','')
        plt.title(hpvName)

        if covMapYmax:
            cov.set_ylim(top=covMapYmax)

        #Plot gene annotations
        glines = []
        glabels = []
        y1end=0
        y2end=0
        ypos1 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/12
        ypos2 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/7.9
        ypos3 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/5.8
        yposlab1 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/8.5
        yposlab2 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/6.2
        yposlab3 = plt.ylim()[0] - (plt.ylim()[1]-plt.ylim()[0])/4.8
        if hpvName in hpvGeneDict:
            ic = 0
            gNameLast = ''
            for gene in hpvGeneDict[hpvName]:
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
                    #y1end = gEnd
                elif gStart >= y2end:
                    ypos = ypos2
                    yposlab = yposlab2
                    #y2end = gEnd
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

        #fig.tight_layout()
        fig.savefig(outputName+'.'+hpvName+'.cov.pdf',bbox_inches='tight',bbox_extra_artists=glines+glabels)
        plt.close(fig)
    return outTable


def main(argv):
    #outTable = mapReads(argv[1:3])
    #print(str([argv[1]])+' '+str(argv[2]))
    #hpvBams, defaultHpvRef=True, hpvRefPath='', filterLowComplex=True, outputName='hpvType', covMapYmax=0
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
