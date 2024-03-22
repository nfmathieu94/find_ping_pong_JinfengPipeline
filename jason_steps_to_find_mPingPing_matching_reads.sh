module load minimap2
module load samtools
minimap2 -a  -x sr left_flank_mPing.fa heg4_fastq/HEG4_2.1_p1.fq.gz -t 24 -o heg4_onlyPing_r1.sam
minimap2 -a  -x sr left_flank_mPing.fa heg4_fastq/HEG4_2.1_p2.fq.gz -t 24 -o heg4_onlyPing_r2.sam

# see reads that are not 101 matches (eg have a flank or a mismatch)
# the -F 4 drops reads that don't align to the genome db at all
samtools view -F 4 heg4_onlyPing_r1.sam | grep -v 101M

# make a BAM file
samtools view -F 4 heg4_onlyPing_r1.sam -OBAM -o match_only.bam
samtools sort match_only.bam match_only.sort.bam

#to visualize
samtools tview match_sort.bam --reference left_flank_mPing.fa

# to mpileup in the context of the reference genome - look for position 16
samtools mpileup match_sort.bam --reference left_flank_mPing.fa
