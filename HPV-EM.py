#!/usr/bin/env python
import argparse as argp
import subprocess as subp
import os
import sys
from whichcraft import which
from EMstep import EmAlgo
from CreateMappedReadTable import mapReads

def prereqs():
    programs = ["python", "bwa", "samtools", "STAR"]
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
            #print("Subprocesss error with code: " +  str(e.returncode))
            exit()
        except:
            print("An unknown error occurred")
            exit()

        sys.stdout = temp

    else:
        try:
            print(' '.join(args))
            subp.check_call(args)
        except subp.CalledProcessError, e:
            #print("Subprocesss error with code: " + str(e.returncode))
            exit()
        except:
            print("An unknown error occurred")
            exit()

    return


def aligntoGenome(args, i, hpvBams):
    cmd(["bwa", 
         "aln", 
         "-t {}".format(args.threads),
         "{}".format(args.reference),
         "{sampleName}.Unmapped.out.mate{i}".format(sampleName=args.outname, i=i)], 
        True, "{sampleName}.{i}.sai".format(sampleName=args.outname, i=i))

    cmd(["bwa", 
         "samse",
         "-n",
         "100",
         "{}".format(args.reference), 
         "{sampleName}.{i}.sai".format(sampleName=args.outname, i=i),
         "{sampleName}.Unmapped.out.mate{i}".format(sampleName=args.outname, i=i)],
        True,"{sampleName}.{i}.aln-se.sam".format(sampleName=args.outname,i=i))

    cmdStr = 'samtools view -F4 -@ {t} -o {sampleName}.aligned.{i}.bam {sampleName}.{i}.aln-se.sam'.format(t=args.threads-1, sampleName=args.outname, i=i)
    cmd(cmdStr.split())
    
    hpvBams.append("{sampleName}.aligned.{i}.bam".format(sampleName=args.outname, i=i))

    if not args.keepint:
        #clean up files
        os.remove("{sampleName}.Unmapped.out.mate{i}".format(sampleName=args.outname, i=i))
        os.remove("{sampleName}.{i}.sai".format(sampleName=args.outname, i=i))
        os.remove("{sampleName}.{i}.aln-se.sam".format(sampleName=args.outname, i=i))

        
def main(): 
    myparse = argp.ArgumentParser(description='Run the HPV-EM genotyping tool', formatter_class=lambda prog: argp.RawTextHelpFormatter(prog, width=99999))
    
    # positional arguments
    myparse.add_argument("reads1", help="single-end FASTQ file or first paired-end FASTQ file")
    myparse.add_argument("reads2", nargs='?', help="(optional) second paired-end FASTQ file", default="not supplied")

    # options
    myparse.add_argument('-t','--threads', type=int,  help="number of threads to use [1]", default=1)
    myparse.add_argument('-r','--reference', help="viral reference genome in FASTA format,\nto be used in place of default HPV reference",default=0)
    myparse.add_argument('-o', '--outname', type=str, help="output file name prefix [./hpvEM]", default='./hpvEM')
    myparse.add_argument('-d', '--disabledust', action='store_true', help="disable filtering of low-complexity reads")
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
    else:
        defaultHpvRef = False

    if not os.path.isfile(args.reference+'.bwt'):
        print 'Indexing viral reference'
        cmd(['bwa',
             'index',
             args.reference])

    if args.threads<1:
        args.threads=1

    if(prereqs() == False):
        sys.exit(1)

    outPath=args.outname.split('/')
    if len(outPath)>1:
        outPath = '/'.join(outPath[:-1])
        if not os.path.isdir(outPath):
            cmd(["mkdir", outPath])

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

        print "Aligning reads to HPV genomes"
        aligntoGenome(args,1,hpvBams)
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

        print "Aligning reads to HPV genomes"
        aligntoGenome(args,1,hpvBams)
        aligntoGenome(args,2,hpvBams)

    if not args.keepint:
        os.remove('{}.Log.progress.out'.format(args.outname))
        os.remove('{}.Log.final.out'.format(args.outname))
        os.remove('{}.Log.out'.format(args.outname))
        os.remove('{}.SJ.out.tab'.format(args.outname))
        os.remove('{}.Chimeric.out.junction'.format(args.outname))
        os.remove('{}.Chimeric.out.sam'.format(args.outname))
        os.remove('{}.Aligned.out.bam'.format(args.outname))

    print "Creating read table"
    readsTable = mapReads(hpvBams, defaultHpvRef=defaultHpvRef, hpvRefPath=args.reference, filterLowComplex=not(args.disabledust), outputName=args.outname)

    print "Running EM algorithm"
    EmAlgo(readsTable, outputName=args.outname, printResult=args.printem)

    sys.exit(0)


if __name__ == '__main__':
    main()
