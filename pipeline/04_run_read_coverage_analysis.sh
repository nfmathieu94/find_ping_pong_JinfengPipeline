#!/usr/bin/bash -l

module load relcate2

python Rice3k_copy_number_depth_window_mPing.py --input Rice3k_3000_RelocaTEi_mPing_NM2 --depth ping_coverage_3k.bam.summary
python Rice3k_copy_number_depth_window_Ping.py --input Rice3k_3000_RelocaTEi_Ping_NM2 --depth ping_coverage_3k.bam.summary
python Rice3k_copy_number_depth_window_Ping.py --input Rice3k_3000_RelocaTEi_Pong_NM2 --depth ping_coverage_3k.bam.summary
