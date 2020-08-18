#!/usr/bin/env python
import argparse as argp
import subprocess as subp
import os
import sys
import shutil
from whichcraft import which
from EMstep import EmAlgo
from CreateMappedReadTable import mapReads

def prereqs():
    programs = ["python", "samtools", "STAR"]
    ready = True;

    for i in range(0, len(programs)):
        if which(programs[i]) is None:
            print(programs[i] + " not installed. Please install " + programs[i])
            ready = False;
    return ready


def cmd(args, write=False, filepath=None):
    if(write==True):
        temp = sys.stdout
        sys.stdout = open(filepath, 'w')

        try:
            subp.check_call(args, stdout=sys.stdout)
        except subp.CalledProcessError, e:
            print("Subprocesss error with code: " +  str(e.returncode))
            sys.exit(e.returncode)
        except:
            print("An unknown error occurred")
            sys.exit(1)

        sys.stdout = temp

    else:
        try:
            print(' '.join(args))
            subp.check_call(args)
        except subp.CalledProcessError, e:
            print("Subprocesss error with code: " + str(e.returncode))
            sys.exit(e.returncode)
        except:
            print("An unknown error occurred")
            exit(1)

    return

        
def main(): 
    myparse = argp.ArgumentParser(description='Run the HPV-EM genotyping tool', formatter_class=lambda prog: argp.RawTextHelpFormatter(prog, width=99999))
    
    # positional arguments
    myparse.add_argument("reads1", help="single-end FASTQ file or first paired-end FASTQ file")
    myparse.add_argument("reads2", nargs='?', help="(optional) second paired-end FASTQ file", default="not supplied")

    # options
    myparse.add_argument('-t','--threads', type=int,  help="number of threads to use [1]", default=1)
    myparse.add_argument('-r','--reference', help="viral reference genome in FASTA format,\nto be used in place of default HPV reference",default=0)
    myparse.add_argument('--starviral', help="path to a directory containing STAR-generated\nviral genome indexes based on the above FASTA",default=0)
    myparse.add_argument('-o', '--outname', type=str, help="output file name prefix [./hpvEM]", default='./hpvEM')
    myparse.add_argument('-d', '--disabledust', action='store_true', help="disable filtering of low-complexity reads")
    myparse.add_argument('--tpm', type=float, help="TPM threshold for identifying a true positive [1.48]", default=1.48)
    myparse.add_argument('-p', '--printem', action='store_true', help="print EM results to STDOUT")
    myparse.add_argument('-k', '--keepint', action='store_true', help="keep intermediate files")

    # other required arguments
    requiredNamed = myparse.add_argument_group('required arguments')
    requiredNamed.add_argument('-s','--stargenome', help="path to a directory containing STAR-generated\nhuman genome indexes", required=True)

    args = myparse.parse_args()

    # finding path to reference
    if(args.reference == 0):
        defaultHpvRef = True
        installDir = os.path.dirname(os.path.abspath(__file__))
        args.reference = installDir+'/reference/combined_pave_hpv.fa'
        args.starviral = installDir+'/reference/combined_pave_hpv_STAR'
    else:
        defaultHpvRef = False
        if args.starviral == 0:
            print('Please provide the path to a folder of STAR indices based on your specified viral genome using the --starviral argument')
            sys.exit(1)

    if args.threads<1:
        args.threads=1

    if(prereqs() == False):
        sys.exit(1)

    outPath=args.outname.split('/')
    if len(outPath)>1:
        outPath = '/'.join(outPath[:-1])
        if not os.path.isdir(outPath):
            cmd(["mkdir", outPath])

    allReadsNum = -1
    hpvBams = []
    if args.reads2 == "not supplied":
        print "Aligning reads to human genome"
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.stargenome),
             "--readFilesIn {}".format(args.reads1),
             "--runThreadN {}".format(args.threads),
             "--chimSegmentMin 18",
             "--outSAMtype BAM Unsorted",
             "--outReadsUnmapped Fastx",
             "--outFilterMultimapNmax 100",
             "--outFileNamePrefix {}.".format(args.outname)])

        with open('{}.Log.final.out'.format(args.outname),'r') as logFile:
            for line in logFile:
                line = line.strip()
                if line.startswith('Number of input reads'):
                    allReadsNum = int(line.split()[-1])
                    print "Total reads: {}".format(allReadsNum)
                    break

        print "Aligning reads to HPV genomes"
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.starviral),
             "--readFilesIn {sampleName}.Unmapped.out.mate1".format(sampleName=args.outname),
             "--runThreadN {}".format(args.threads),
             "--twopassMode Basic",
             "--outSAMtype BAM Unsorted",
             "--outSAMattributes NH HI NM MD AS XS",
             "--outFilterMultimapNmax 999",
             "--outFilterMismatchNmax 999",
             "--outFilterMismatchNoverLmax 0.08",
             "--outFileNamePrefix {}.1.".format(args.outname)])
        hpvBams.append('{}.1.Aligned.out.bam'.format(args.outname))

    else:
        print "Aligning reads to human genome"
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.stargenome),
             "--readFilesIn {} {}".format(args.reads1, args.reads2),
             "--runThreadN {}".format(args.threads),
             "--chimSegmentMin 18",
             "--outSAMtype BAM Unsorted",
             "--outReadsUnmapped Fastx",
             "--outFilterMultimapNmax 100",
             "--outFileNamePrefix {}.".format(args.outname)])

        with open('{}.Log.final.out'.format(args.outname),'r') as logFile:
            for line in logFile:
                line = line.strip()
                if line.startswith('Number of input reads'):
                    allReadsNum = int(line.split()[-1])
                    break

        print "Aligning reads to HPV genomes"
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.starviral),
             "--readFilesIn {sampleName}.Unmapped.out.mate1".format(sampleName=args.outname),
             "--runThreadN {}".format(args.threads),
             "--twopassMode Basic",
             "--outSAMtype BAM Unsorted",
             "--outSAMattributes NH HI NM MD AS XS",
             "--outFilterMultimapNmax 999",
             "--outFilterMismatchNmax 999",
             "--outFilterMismatchNoverLmax 0.08",
             "--outFileNamePrefix {}.1.".format(args.outname)])
        hpvBams.append('{}.1.Aligned.out.bam'.format(args.outname))
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.starviral),
             "--readFilesIn {sampleName}.Unmapped.out.mate2".format(sampleName=args.outname),
             "--runThreadN {}".format(args.threads),
             "--twopassMode Basic",
             "--outSAMtype BAM Unsorted",
             "--outSAMattributes NH HI NM MD AS XS",
             "--outFilterMultimapNmax 999",
             "--outFilterMismatchNmax 999",
             "--outFilterMismatchNoverLmax 0.08",
             "--outFileNamePrefix {}.2.".format(args.outname)])
        hpvBams.append('{}.2.Aligned.out.bam'.format(args.outname))

    if not args.keepint:
        os.remove('{}.Log.progress.out'.format(args.outname))
        os.remove('{}.Log.final.out'.format(args.outname))
        os.remove('{}.Log.out'.format(args.outname))
        os.remove('{}.SJ.out.tab'.format(args.outname))
        os.remove('{}.Chimeric.out.junction'.format(args.outname))
        os.remove('{}.Aligned.out.bam'.format(args.outname))
        os.remove('{}.Unmapped.out.mate1'.format(args.outname))

        os.remove('{}.1.Log.progress.out'.format(args.outname))
        os.remove('{}.1.Log.final.out'.format(args.outname))
        os.remove('{}.1.Log.out'.format(args.outname))
        os.remove('{}.1.SJ.out.tab'.format(args.outname))
        shutil.rmtree('{}.1._STARgenome'.format(args.outname))
        shutil.rmtree('{}.1._STARpass1'.format(args.outname))
        if args.reads2 != "not supplied":
            os.remove('{}.Unmapped.out.mate2'.format(args.outname))
            os.remove('{}.2.Log.progress.out'.format(args.outname))
            os.remove('{}.2.Log.final.out'.format(args.outname))
            os.remove('{}.2.Log.out'.format(args.outname))
            os.remove('{}.2.SJ.out.tab'.format(args.outname))
            shutil.rmtree('{}.2._STARgenome'.format(args.outname))
            shutil.rmtree('{}.2._STARpass1'.format(args.outname))

    print "Creating read table"
    readsTable = mapReads(hpvBams, defaultHpvRef=defaultHpvRef, hpvRefPath=args.reference, filterLowComplex=not(args.disabledust), outputName=args.outname)

    print "Running EM algorithm"
    EmAlgo(readsTable, allReadsNum, thresholdTpm=args.tpm, outputName=args.outname, printResult=args.printem)

    sys.exit(0)


if __name__ == '__main__':
    main()
