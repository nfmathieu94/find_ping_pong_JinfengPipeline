# Determining Ping Locations

Information from email to Jinfeng:  

 
Sorry for the delay. Identifying ping and pong is not an easy job. It depends on the situations.  
In the RILs, because we know the sequences of Ping and mPing in HEG4 and NB. It is much easier.  
But we still use my pipeline to identify the potential Ping and Pong sites.  
And Lu Lu did PCR to confirm new Ping sites that are not present in the parental strains (HEG4 and NB).  
In rice 3000 accessions, because the Ping sequences in each accessions might be quite different.  
So I run my pipeline to identify Ping and Pong sites. I also run analysis based on read depth to estimate  
Ping and Pong copy numbers to make sure the results are robust (Supplementary Figure 1 in Nature communication paper).  
We did not do PCR confirmation in rice 3000 accessions project as we do not have materials. 
 
The pipeline I mentioned is described [here](https://github.com/stajichlab/Dynamic_rice_publications/blob/master/rice_3k_mPing_scripts/Copy_numbers_characterization/work.sh)

## Protocol

The following steps are from the suggested protocol in Jinfeng's email.

### Step1

Step1, Run RelocaTE2 to detect mPing, Ping, Pong insertion, respectively. The results are in three different folder.

python ReNameSRA_RelocaTEi_mPing.py --input fastq --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa > log 2>&1  

python ReNameSRA_RelocaTEi_mPing.py --input fastq --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/ping.fa > log 2>&1  

python ReNameSRA_RelocaTEi_mPing.py --input fastq --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/pong.fa > log 2>&1  

### Step2

Step2, Run scripts to identify Ping, Pong, and mPing sequentially based on RelocaTE2 results folders generated above.

python ReNameSRA_sum_Ping.py --input fastq_Ping  

python ReNameSRA_sum_Pong.py --input fastq_Pong  

python ReNameSRA_sum_mPing.py --input fastq_mPing  

### Step3

Step3, Map mPing/Ping/Pong reads back to mPing/Ping/Pong elements and summary relative depth of mPing/Ping/Pong to actin or some other CONTROL regions.  

This will give you an estimate of mPing/Ping/Pong copy numbers. I did not find the commend how I generate the mpileup file.  
But the idea is to use RelocaTE2 clipped reads that only matched to mPing/Ping/Pong elements to mapped back to mPing/Ping/Pong elements.  
This gives you accurate alignments than mapping raw reads. 

python Rice3k_copy_number_depth_window_mPing.py --input Rice3k_3000_RelocaTEi_mPing_NM2 --depth ping_coverage_3k.bam.summary

python Rice3k_copy_number_depth_window_Ping.py --input Rice3k_3000_RelocaTEi_Ping_NM2 --depth ping_coverage_3k.bam.summary

python Rice3k_copy_number_depth_window_Ping.py --input Rice3k_3000_RelocaTEi_Pong_NM2 --depth ping_coverage_3k.bam.summary

## Testing

Currently this procedure is being tested by using HEG4 short read data. HEG4 was chosen because we know where mPing insertions
are, so this can serve as a true positive.

### Pipeline Directory

This directory has the shell scripts that run the different python code for each of the steps above.

Some of these scripts are not fully filled out because issues are occuring with the first step.


### Current Issues

An issue is occuring when running ./pipeline/01_run_relocate_mping_ping_pong.sh

The scripts generate output directories that are empty, and the log file output is not informative.


The scripts suggested to run have hard coded paths that do not exist anymore.  
The scripts also seem to use an old queing system, so these need to be updated to slurm.



