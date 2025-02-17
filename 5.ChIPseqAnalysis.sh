# ChIP-seq Analysis

fq1=ND1_ChIP_1.fq.gz
fq2=ND1_ChIP_2.fq.gz
sam=ND1.sam
bam=ND1.bam

# QC
fastqc -o ND1 -f fastq $fq1
fastqc -o ND1 -f fastq $fq2
trim_galore -q 20 --phred33 --paired $fq1 $fq2 --gzip -o ND1 

# Alignment
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -p 10 -q \
	-x RatBN7.2 \
	-1 $fq1-2 $fq2 \
	-S $sam

samtools view -Sb -@ 4 $sam | samtools sort -@ 4 -m 4g > $bam

# Remove Duplicates
picard MarkDuplicates \
	REMOVE_DUPLICATES=true \
	I=$bam O=ND1_dedip.bam 

# call peak
macs2 callpeak -t ND1_dedip.bam  -f BAMPE -q 0.05 -B --keep-dup 1 --outdir ND1

# Peak visulization using IGV

# peak annotation
annotatePeaks.pl ND1_peaks.narrowPeak Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa -gtf RatBN7.2.gtf -gid > ND1_peaks.narrowPeak.homeranno.txt

# signal distribution
BLACKLIST="rn7_liftOver_mm10-blacklist.v2.nochr.bed"
bed=NeuroD1_merged_peak_summits.bed
input_bw=(D0_H3K27ac_treat_pileup.bw  D1_H3K27ac_treat_pileup.bw  D2_H3K27ac_treat_pileup.bw  D5_H3K27ac_treat_pileup.bw)
labels=(D0_H3K27ac  D1_H3K27ac  D2_H3K27ac  D5_H3K27ac)

computeMatrix reference-point -p 10 -R ${bed} -S ${input_bw} --samplesLabel ${labels} \
    --referencePoint center --upstream 2000 --downstream 2000 --missingDataAsZero --skipZeros \
    --blackListFileName ${BLACKLIST} --outFileName center.mat.gz\
    --outFileNameMatrix center.values.tab --outFileSortedRegions center.sorted.bed

plotProfile -m center.mat.gz -out center.pdf 

