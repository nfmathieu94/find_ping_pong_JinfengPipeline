#!/usr/bin/env python2
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import time
from Bio import SeqIO

script_dir = os.path.dirname(os.path.realpath(__file__))

def usage():
    test="name"
    message='''
python ReNameSRA_RelocaTEi.py --input Japonica_fastq

Run RelocaTEi for rice strain in Japonica_fastq

    '''
    print(message)


def runjob(script, lines):
    cmd = 'sbatch --array=1-%s --nodes=1 --ntasks-per-node=1 --time=100:00:00 --mem=10G %s' % (lines, script)
    #print cmd 
    os.system(cmd)


def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = '%s_RelocaTEi' %(os.path.abspath(args.input))

    if not args.genome:
        #args.genome = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa'
        args.genome = os.path.join('genome', 'MSU_r7.fa')
  
    if not args.repeat:
        #args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa'
        args.repeat = os.path.join('lib', 'mping.fa')


    #RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2_mPing/scripts/relocaTE.py'
    RelocaTE = os.path.join('relocate', 'relocaTE.py')



    Reference= os.path.abspath(args.genome)
    Repeat   = os.path.abspath(args.repeat)
    project = os.path.split(args.output)[1]
    cpu = 16
    
    if not os.path.exists(project):
        os.mkdir(project)
    print(project)

    read_dirs = glob.glob('%s/*' %(os.path.abspath(args.input)))
    ofile = open('%s.run.sh' %(args.output), 'w')

    for read_dir in sorted(read_dirs):
        outdir = '%s/%s_RelocaTEi' %(os.path.abspath(args.output), os.path.split(read_dir)[1])
        existingTE  = '%s.mPing.RepeatMasker.out' %(Reference)
        # relocate will not run if there is result exists
        if not os.path.exists(outdir):
        #if 1:
            relocaTE = '%s --te_fasta %s --genome_fasta %s --fq_dir %s --outdir %s --reference_ins %s' %(RelocaTE, Repeat, Reference, read_dir, outdir, existingTE)
            shell    = 'bash %s/run_these_jobs.sh > %s/run.log 2> %s/run.log2' %(outdir, outdir, outdir)
            os.system(relocaTE)
            #print >> ofile, relocaTE
            print >> ofile, shell
    ofile.close()
    runjob('%s.run.sh' %(args.output), 1)

    for read_dir in sorted(read_dirs):
        outdir = '%s/%s_RelocaTEi' % (os.path.abspath(args.output), os.path.split(read_dir)[1])
        existingTE = '%s.mPing.RepeatMasker.out' % Reference
    if not os.path.exists(outdir):
        print("Creating output directory:", outdir)  # Debugging print statement
        relocaTE = '%s --te_fasta %s --genome_fasta %s --fq_dir %s --outdir %s --reference_ins %s' % (
            RelocaTE, Repeat, Reference, read_dir, outdir, existingTE)
        shell = 'bash %s/run_these_jobs.sh > %s/run.log 2> %s/run.log2' % (outdir, outdir, outdir)
        os.system(relocaTE)
        print("Shell command:", shell)  # Debugging print statement
        print >> ofile, shell  # Writing to the batch script file

    
    ofile.close()
    runjob('%s.run.sh' %(args.output), 1)
 
if __name__ == '__main__':
    main()

