#!/usr/bin/bash -l

module load relocate2

# Path to step1 script
reloc_script=dynamic_rice_find_ping_scripts/testing.py
input=element_sim_fastq
genome=genome/Chr4.fasta
mping=lib/mping.fa
ping=lib/elements/ping.fa
pong=lib/elements/pong.fa


# Only running for mPing until things are working
python2 $reloc_script  --input $input  --repeat $mping --genome $genome

#python2 $RELOC_SCRIPT -input fastq --repeat ./lib/elements/ping.fa

#python2 $RELOC_SCRIPT -input fastq --repeat ./lib/elements/pong.fa



