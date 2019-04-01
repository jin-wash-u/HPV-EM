import argparse as argp
import subprocess as subp
import os
import sys
import glob
from whichcraft import which


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
            subp.check_call(args)
        except subp.CalledProcessError, e:
            #print("Subprocesss error with code: " + str(e.returncode))
            exit()
        except:
            print("An unknown error occurred")
            exit()

    return



def aligntoGenome(nameOnly, i, args, topdirectory):
    cmd(
        ["bwa", 
        "aln", 
        "-t {}".format(args.threads),
        "{}".format(args.reference),
        "{sampleName}.Unmapped.out.mate{i}".format(sampleName=nameOnly, i=i)], 
        True, "{sampleName}.{i}.sai".format(sampleName=nameOnly, i=i))

    cmd(
        ["bwa", 
        "samse", 
        "{}".format(args.reference), 
        "{sampleName}.{i}.sai".format(sampleName=nameOnly, i=i),
        "{sampleName}.Unmapped.out.mate{i}".format(sampleName=nameOnly, i=i)],
        True,"{sampleName}.{i}.aln-se.sam".format(sampleName=nameOnly,i=i))

    cmd(
        ["samtools",
        "view",
        "-F4"
        "-@ {}".format((args.threads)),
        "{sampleName}.{i}.aln-se.sam".format(sampleName=nameOnly, i=i)], 
        True,"{topdir}/{sampleName}/HPV.aligned.{i}.sam".format(topdir=topdirectory, sampleName=nameOnly, i=i))

    for format in {'*.sam','*.sai'}:
        for file in glob.glob(format):
            os.remove(file)



def main(): 

    myparse = argp.ArgumentParser(description='Runs the HPV alignment tool')
    
    # positional arguments
    myparse.add_argument("reads1", help="single end FASTQ file or first paired end FASTQ file")
    myparse.add_argument("reads2", nargs='?', help="(optional) second paired end FASTQ file", default="not supplied")


    # other required arguments
    myparse.add_argument('-s','--stargenome', required=True, type=dir, help="path to a directory containing STAR generated human genome files")


    # options
    myparse.add_argument('-t','--threads', type=int,  help="number of threads to use [1]", default=1)
    myparse.add_argument('-r','--reference', help="viral reference genome in FASTA format, to be used in place of default HPV reference",default=0)
    myparse.add_argument('-o', '--outname', type=str, help="output file name prefix [./hpvEM]", default='./hpvEM')
    myparse.add_argument('-d', '--dust', type=int, help="0 to disable filtering of low-complexity reads", default=1)
    


    args = myparse.parse_args()

    # finding path to reference
    if(args.reference == 0):
        defaultHpvRef = True
    else:
        defaultHpvRef = False
        path = args.reference


    if(prereqs() == False):
        exit()

    topdirectory = os.getcwd()

    nameOnly = ('.').join(args.reads1.split('.')[:-1]) # extracting name of file without extension

    if (not(os.path.isdir("{sampleName}".format(sampleName=nameOnly)))):
        cmd(["mkdir", nameOnly])
    else:
        num = 1

        while(os.path.isdir("{sampleName}_{number}".format(sampleName=nameOnly, number=num))):
            num = num + 1

        cmd(["mkdir", "{sampleName}_{number}".format(sampleName=nameOnly, number=num)])


    if((args.reads1.lower().endswith(".bam")) or (args.reads2 != "not supplied")):
        
        if(args.reads1.lower().endswith(".bam")): # if bam file given as input, convert to fastq files
            
            cmd(["echo", "Extracting raw reads"]) 
            cmd(["samtools", "fastq",
                 "-1{}.1.fq".format(nameOnly),
                 "-2{}.2.fq".format(nameOnly),
                 "-0{}".format(os.devnull), 
                 "-n", "-F 0x900", "-@ {}".format(args.threads),
                 "{}".format(args.reads1)])

            firstSample = nameOnly + ".1"
            secondSample = nameOnly + ".2"

        else:
            firstSample = nameOnly
            secondSample = ('.').join(args.reads2.split('.')[:-1])
        
        cmd(["echo", "Aligning reads to human genome"])


        cmd(
            ["STAR", 
            "--genomeDir {path}".format(path=args.stargenome),
            "--readFilesIn {firstSample}.fq {secondSample}.fq".format(firstSample=firstSample, secondSample=secondSample),
            "--runThreadN {}".format(args.threads),
            "--chimSegmentMin 18",
            "--outSAMtype BAM Unsorted",
            "--outReadsUnmapped Fastx",
            "--outFilterMultimapNmax 100",
            "--outFileNamePrefix ./{}.".format(nameOnly)])

        if(args.reads1.lower().endswith(".bam")):
            cmd(
            ["rm", 
            "{}.1.fq".format(nameOnly), 
            "{}.2.fq".format(nameOnly)])


        cmd(["echo", "Aligning reads to HPV genomes"])

        aligntoGenome(nameOnly,1,args,topdirectory)
        aligntoGenome(nameOnly,2,args,topdirectory)


    else:
        cmd(["echo", "Aligning reads to human genome"])


        cmd(
            ["STAR", 
            "--genomeDir {path}".format(path=args.stargenome),
            "--readFilesIn {sampleName}.fq".format(sampleName=nameOnly),
            "--runThreadN {}".format(args.threads),
            "--chimSegmentMin 18",
            "--outSAMtype BAM Unsorted",
            "--outReadsUnmapped Fastx",
            "--outFilterMultimapNmax 100",
            "--outFileNamePrefix ./{}.".format(args.outname)])

        cmd(["echo", "Aligning reads to HPV genomes"])

        aligntoGenome(nameOnly,1,args,topdirectory)



    cmd(["echo", "Required file clean up"])

    for format in {'*.out','*.junction','*.tab','*.Aligned.*', '*.mate*'}:
        for file in glob.glob(format):
            os.remove(file)

    exit()


if __name__ == '__main__':
    main()
