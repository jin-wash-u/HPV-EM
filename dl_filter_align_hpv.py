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
        "-t {}".format(args.cpus),
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
        "-@ {}".format((args.cpus - 1)),
        "{sampleName}.{i}.aln-se.sam".format(sampleName=nameOnly, i=i)], 
        True,"{topdir}/{sampleName}/HPV.aligned.{i}.sam".format(topdir=topdirectory, sampleName=nameOnly, i=i))

    for format in {'*.sam','*.sai'}:
        for file in glob.glob(format):
            os.remove(file)



def main(): 

    myparse = argp.ArgumentParser(description='Runs the HPV alignment tool')
    myparse.add_argument('sampleName', metavar="sampleName", help="name of the first sample to be aligned")
    myparse.add_argument('reference', metavar="refFasta", help="the reference .fa file for the program")
    myparse.add_argument('path', metavar="path", help="path to Human Genome (.fa, .fasta) file")
    myparse.add_argument('-@', dest="cpus", type=int, default=2, help="number of CPUS for processing")
    myparse.add_argument('-2', dest="otherSample", metavar="otherSample", default="not supplied", help="name of the (optional) second sample to be aligned")


    args = myparse.parse_args()

    if(prereqs() == False):
        exit()

    topdirectory = os.getcwd()

    nameOnly = ('.').join(args.sampleName.split('.')[:-1])

    if (not(os.path.isdir("{sampleName}".format(sampleName=nameOnly)))):
        cmd(["mkdir", nameOnly])
    else:
        num = 1

        while(os.path.isdir("{sampleName}_{number}".format(sampleName=nameOnly, number=num))):
            num = num + 1

        cmd(["mkdir", "{sampleName}_{number}".format(sampleName=nameOnly, number=num)])


    if((args.sampleName.lower().endswith(".bam")) or (args.otherSample != "not supplied")):
        
        if(args.sampleName.lower().endswith(".bam")): # if bam file given as input, convert to fastq files
            
            cmd(["echo", "Extracting raw reads"]) 
            cmd(["samtools", "fastq",
                 "-1{}.1.fq".format(nameOnly),
                 "-2{}.2.fq".format(nameOnly),
                 "-0{}".format(os.devnull), 
                 "-n", "-F 0x900", "-@ {}".format(args.cpus-1),
                 "{}".format(args.sampleName)])

            firstSample = nameOnly + ".1"
            secondSample = nameOnly + ".2"

        else:
            firstSample = nameOnly
            secondSample = ('.').join(args.otherSample.split('.')[:-1])
        
        cmd(["echo", "Aligning reads to human genome"])


        cmd(
            ["STAR", 
            "--genomeDir {path}".format(path=args.path),
            "--readFilesIn {firstSample}.fq {secondSample}.fq".format(firstSample=firstSample, secondSample=secondSample),
            "--runThreadN {}".format(args.cpus),
            "--chimSegmentMin 18",
            "--outSAMtype BAM Unsorted",
            "--outReadsUnmapped Fastx",
            "--outFilterMultimapNmax 100",
            "--outFileNamePrefix ./{}.".format(nameOnly)])

        if(args.sampleName.lower().endswith(".bam")):
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
            "--genomeDir {path}".format(path=args.path),
            "--readFilesIn {sampleName}.fq".format(sampleName=nameOnly),
            "--runThreadN {}".format(args.cpus),
            "--chimSegmentMin 18",
            "--outSAMtype BAM Unsorted",
            "--outReadsUnmapped Fastx",
            "--outFilterMultimapNmax 100",
            "--outFileNamePrefix ./{}.".format(nameOnly)])

        cmd(["echo", "Aligning reads to HPV genomes"])

        aligntoGenome(nameOnly,1,args,topdirectory)



    cmd(["echo", "Required file clean up"])

    for format in {'*.out','*.junction','*.tab','*.Aligned.*', '*.mate*'}:
        for file in glob.glob(format):
            os.remove(file)

    exit()


if __name__ == '__main__':
    main()
