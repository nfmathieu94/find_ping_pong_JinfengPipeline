#!/usr/bin/bash -l

module load relcate2

# Path to step1 script
RELOC_SCRIPT=./dynamic_rice_find_ping_scripts/ReNameSRA_RelocaTEi_mPing.py


# Only running for mPing until things are working
python2 $RELOC_SCRIPT  -input heg4_fastq  --repeat ./lib/elements/mping.fa

#python2 $RELOC_SCRIPT -input fastq --repeat ./lib/elements/ping.fa

#python2 $RELOC_SCRIPT -input fastq --repeat ./lib/elements/pong.fa



