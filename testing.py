#!/usr/bin/env python2
import sys
from collections import defaultdict
import re
import os
import argparse
import glob
from Bio import SeqIO


def usage():
    test = "name"
    message = '''
python ReNameSRA_RelocaTEi.py --input Japonica_fastq

Run RelocaTEi for rice strain in Japonica_fastq
    '''
    print(message)


def runjob(script, lines):
    cmd = 'sbatch --array=1-%s --nodes=1 --ntasks-per-node=1 --time=100:00:00 --mem=10G %s' % (lines, script)
    os.system(cmd)


def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile, "fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open(infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line)
                if unit[0] not in data:
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
        assert args.input, "Input directory is required"
    except AssertionError as e:
        print(e)
        usage()
        sys.exit(2)

    if not args.output:
        args.output = os.path.abspath('%s_RelocaTEi' % args.input)

    if not args.genome:
        args.genome = os.path.join('genome', 'MSU_r7.fa')

    if not args.repeat:
        args.repeat = os.path.join('lib', 'mping.fa')

    RelocaTE = os.path.join('relocate', 'relocaTE.py')
    Reference = os.path.abspath(args.genome)
    Repeat = os.path.abspath(args.repeat)
    project = os.path.split(args.output)[1]

    if not os.path.exists(project):
        os.mkdir(project)

    read_dirs = glob.glob('%s/*' % os.path.abspath(args.input))

    with open('%s.run.sh' % args.output, 'w') as ofile:
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
                print >> ofile, "Logging:", outdir
                print >> ofile, "Logging:", shell


if __name__ == '__main__':
    main()
