#!/usr/bin/bash -l

module load relocate2

# Will modify this scripts once ./00_run_relocate_mping_ping_pong.sh works


python2 ReNameSRA_sum_Ping.py --input fastq_Ping
python2 ReNameSRA_sum_Pong.py --input fastq_Pong
python2 ReNameSRA_sum_mPing.py --input fastq_mPing
