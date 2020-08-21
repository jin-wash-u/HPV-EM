#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import math
import sys
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Difference value between successive loglikelihoods below which algorithm is deemed to have converged
conVal = 1e-4
numIter = 5
maxSteps = 500

class mappedRead:
    def __init__(self, inputList):
        self.hpvType = inputList[0]
        self.readInfo = []
        self.readNum = 0
        self.genes = []
        count=0
        for val in inputList[1:]:
            if (count%4)==0:
                val = int(val)
                self.readInfo.append([val])
                self.readNum+=val
            elif (count%4)==3:
                if not val:
                    self.genes.append([])
                else:
                    self.genes.append(val.split(','))
            else:
                self.readInfo[-1].append(float(val))
            count+=1

            
def natural_order(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key[1]) ]
    return [i[0] for i in sorted(enumerate(l), key = alphanum_key)]


def EmAlgo(readsTable, allReadsNum, thresholdTpm=1.5, outputName='hpvType', printResult=True):
    mappedReads = []
    totalReads = 0
    uniqReads = int(readsTable[1].split('\t')[0])
    isReadAmbig = readsTable[1].split('\t')[1:]
    for line in readsTable[2:]:
        line = line.split('\t')
        mappedReads.append(mappedRead(line))
        totalReads+=mappedReads[-1].readNum

    if mappedReads:
        # Initialize EM algorithm
        m = len(mappedReads[0].readInfo) #number of reads
        k = len(mappedReads) #number of HPV types

        #Parameters
        for ni in range(numIter):
            if ni==0:
                lOut = -float('inf')
                err = 0.005
                phi = [1./k]*k
            else:
                err = 0.05*np.random.random(1)[0]
                phi = np.random.random(k)
                phi /= phi.sum()

            w=np.zeros([m,k])
            steps=0

            # Calculate initial l
            l=-float('inf')

            converged=False
            while not converged:
                steps+=1
                if(steps>maxSteps):
                    raise RuntimeError('EM algorithm failed to converge after {} steps; aborting.'.format(maxSteps))
                    
                # E step
                for j,hpv in enumerate(mappedReads):
                    for i,readInfo in enumerate(hpv.readInfo):
                        if readInfo[0]:
                            Lm = readInfo[1]
                            Le = readInfo[2]
                            w[i,j] = ((1.-err)**Lm * err**Le) * phi[j]
                for i in range(m):
                    w[i,:] = w[i,:]/sum(w[i,:])

                # M step
                ## err
                Bnum = 0
                Bden = 0
                for j,hpv in enumerate(mappedReads):
                    for i,readInfo in enumerate(hpv.readInfo):
                        if readInfo[0]:
                            Lm = readInfo[1]
                            Le = readInfo[2]
                            Bnum += w[i,j]*Le
                            Bden += w[i,j]*Lm
                B = Bnum/Bden
                err = B/(1.+B)

                ## phi
                for j in range(k):
                    phi[j] = sum(w[:,j])/m

                # Calculate loglikelihood, check change
                l0 = l
                l = 0
                for j,hpv in enumerate(mappedReads):
                    for i,readInfo in enumerate(hpv.readInfo):
                        if readInfo[0]:
                            Lm = readInfo[1]
                            Le = readInfo[2]
                            l +=  w[i,j] * math.log(((1.-err)**Lm * err**Le * phi[j])/w[i,j])

                if (l-l0) < conVal:
                    converged=True
                    if ni==0 or l>lOut:
                        lOut = l
                        errOut = err
                        phiOut = phi
                        stepsOut = steps
                        iterOut = ni

        # Print out results:
        types = []
        typesAll = []
        readProps = []
        emProps = []
        output=[]
        hpvGeneReadCountsDict = {}
        geneNamesSet = set()
        
        # Get number of reads that pass TPM threshold:
        totalOutReads = 0
        for j,hpv in enumerate(mappedReads):
            if uniqReads*phiOut[j]*1e6/allReadsNum > thresholdTpm:
                totalOutReads += int(round(uniqReads*phiOut[j]))
            
        for j,hpv in enumerate(mappedReads):
            hpvName = hpv.hpvType
            if uniqReads*phiOut[j]*1e6/allReadsNum > thresholdTpm:
                types.append(hpvName)
                typesAll.append(hpvName)
                output.append('{!s}\t{:d}\t{:.5f}\t{:d}\t{:.5f}'.format(hpvName,
                                    hpv.readNum, float(hpv.readNum)/totalReads,
                                    int(round(uniqReads*phiOut[j])), round(uniqReads*phiOut[j])/totalOutReads))
                readProps.append(float(hpv.readNum)/totalReads)

                emProps.append(round(uniqReads*phiOut[j])/totalOutReads)

                # Get per-gene read counts
                if hpvName not in hpvGeneReadCountsDict:
                    hpvGeneReadCountsDict[hpvName] = {}

                for ii in range(len(hpv.readInfo)):
                    try:
                        geneList = hpv.genes[ii]
                    except:
                        print("Len hpv.readInfo: "+str(len(hpv.readInfo)))
                        print("current index: "+str(ii))
                        print("error from "+str(hpv.genes))
                        sys.exit(1)
                    if isReadAmbig[ii] == 'U':
                        val = 1
                    else:
                        val = phiOut[j]
                    if geneList:
                        for gene in geneList:
                            geneNamesSet.add(gene)
                            if gene not in hpvGeneReadCountsDict[hpvName]:
                                hpvGeneReadCountsDict[hpvName][gene] = val
                            else:
                                hpvGeneReadCountsDict[hpvName][gene] += val
            else:
                if os.path.exists(outputName+'.'+hpvName+'.cov.pdf'):
                    os.remove(outputName+'.'+hpvName+'.cov.pdf')
                if float(hpv.readNum)/totalReads > 1e-4:
                    typesAll.append(hpvName)
                    readProps.append(float(hpv.readNum)/totalReads)
                    emProps.append(0.0)

               
        # Print and write results to output file
        if output:
            lOrd = natural_order(types)
            lOrdAll = natural_order(typesAll)
            if printResult:
                print('Converged to < {:.1e} in {:d} iterations'.format(conVal, stepsOut))
                print('err\t{:.5f}'.format(errOut))
                print('HPVtype\tMappedReads\tMappedProportion\tMLE_Reads\tMLE_Probability')
                for i in lOrd:
                    print(output[i])

            with open(outputName+'.results.tsv','w') as fOut:
                fOut.write('Converged to < {:.1e} in {:d} iterations\n'.format(conVal, stepsOut))
                fOut.write('err\t{:.5f}\n'.format(errOut))
                fOut.write('HPVtype\tMappedReads\tMappedProportion\tMLE_Reads\tMLE_Probability\n')
                for i in lOrd:
                    fOut.write(output[i]+'\n')           

            # Write out read counts table
            geneNamesList = sorted(geneNamesSet)
            with open(outputName+'.readCounts.tsv','w') as fCounts:
                fCounts.write('Type\t'+'\t'.join(geneNamesList)+'\n')
                for hpvType in hpvGeneReadCountsDict:
                    fCounts.write(hpvType)
                    for gene in geneNamesList:
                        if gene in hpvGeneReadCountsDict[hpvType]:
                            fCounts.write('\t{:.3f}'.format(hpvGeneReadCountsDict[hpvType][gene]))
                        else:
                            fCounts.write('\t0')
                    fCounts.write('\n')

            # Plot pie charts of probabilities, before and after
            ordTypes = map(typesAll.__getitem__, lOrdAll)
            ordReadProps = map(readProps.__getitem__, lOrdAll)
            ordEmProps = map(emProps.__getitem__, lOrdAll)
            ordRpLabels = ['{}\n{:.1f}%'.format(typ, ordReadProps[i]*100) if ordReadProps[i] > 0.01 else '' for i,typ in enumerate(ordTypes)]
            ordEmLabels = ['{}\n{:.1f}%'.format(typ, ordEmProps[i]*100) if ordEmProps[i] > 0.01 else '' for i,typ in enumerate(ordTypes)]

            # Create custom color map for pie charts
            fig, axs = plt.subplots(1, 2, figsize=(12, 6), subplot_kw=dict(aspect="equal"))
            axs[0].set_prop_cycle('color', plt.cm.gist_rainbow(np.linspace(0,1,len(typesAll))))
            axs[1].set_prop_cycle('color', plt.cm.gist_rainbow(np.linspace(0,1,len(typesAll))))
            wedges, texts = axs[0].pie(ordReadProps)
            for i, p in enumerate(wedges):
                ang = (p.theta2 - p.theta1)/2. + p.theta1
                y = 0.6*np.sin(np.deg2rad(ang))
                x = 0.6*np.cos(np.deg2rad(ang))
                axs[0].annotate(ordRpLabels[i], xy=(x, y), ha='center', va='center')

            wedges1, texts1 = axs[1].pie(ordEmProps)
            for i, p in enumerate(wedges1):
                ang = (p.theta2 - p.theta1)/2. + p.theta1
                y = 0.6*np.sin(np.deg2rad(ang))
                x = 0.6*np.cos(np.deg2rad(ang))
                axs[1].annotate(ordEmLabels[i], xy=(x, y), ha='center', va='center')

            axs[1].legend(wedges, ordTypes, title="Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

            axs[0].set_title("Mapped Read Proportions")
            axs[1].set_title("Maximum Likelihood Estimate")

            fig.subplots_adjust(wspace = 0, right = 0.8)
            fig.tight_layout(rect=[0, 0, 0.9, 0.9])

            fig.savefig(outputName+'.props.pdf')
            plt.close(fig)

        else:
            with open(outputName+'.results.tsv','w') as fOut:
                fOut.write('No HPV types detected\n')
            if printResult:
                print('No HPV types detected')
    else:
        with open(outputName+'.results.tsv','w') as fOut:
            fOut.write('No HPV types detected\n')
        if printResult:
            print('No HPV types detected')


def main(argv):
    readTableFile = sys.argv[1]
    readTable = []
    if len(sys.argv) > 3:
        outputName = sys.argv[3]
    else:
        outputName = 'hpvType'
    with open(readTableFile,'r') as inFile:
        for line in inFile:
            readTable.append(line.strip('\n'))
    EmAlgo(readTable, allReadsNum=int(sys.argv[2]), thresholdTpm=1.48, outputName=outputName, printResult=True)

        
if __name__=="__main__":
    main(sys.argv)
