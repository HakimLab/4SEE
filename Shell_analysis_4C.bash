
#alignment using bowtie

bowtie -m 1 -q -S --trim5 NumToRemove Genome_indexs Sample.fastq Sample.sam > & log_Sample &

# NumToRemove - The length of the primer that is in the begging of the sequence that needs to be removed before alignment. The RE site should be remain.
# Genome_indexs - The bowtie indices for the genome
# Sample.fastq - input fastq file with sequences of one bait in one condition.
# Sample.sam - output file in sam format.

#Convert the sam file to bed format for further use.
sam2bed Sample.sam Sample.bed #sam2bed is a script in python.
# A conversion to bam is also possible
samtools view -S -b Sample.sam > Sample.bam

#Assign reads to RE sites
intersectBed -a RE_sites_genome.bed -b Sample.bed -c > RE_assgin_Sample.bedgraph

# RE_sites_genome.bed - input file with all RE sites in the genome in a bed format - chr start end. The length of the intervals should be the length of the RE sequence (4 or 6)
# RE_assgin_Sample.bedgraph - output file in bedgraph format. The file can be use in different genome browsers.

#For sgr file format (that can also be used in IGV genome browser) -
cat RE_assgin_Sample.bedgraph | cut -f1,2,4 > RE_assgin_Sample.sgr

### P value analysis
# The main script is in R language - Pval_analysis_4C.r
# The script required genomic windows spanning the RE sites. The length of the window can change according to the analysis (from 50kb-1Mb).

#Create genomic windows - 

awk -v OFS='\t' '{if (($2-N/2)>0) print $1,$2-N/2,$3+N/2; else print $1,1,$3+N/2 }' RE_sites_genome.bed | \
	intersectBed -a stdin -b genome_file.bed > RE_sites_genome_N_windows.bed


# N - The length of the window size. For example if the length of the window, N, is 100kb, N/2 will be 50kb.
# genome_file.bed - Chromosome coordinates of the genome in bed format - chr 1 length of chromosome.

#Count the total number of RE sites in each window
intersectBed -a RE_sites_genome_N_windows.bed -b RE_sites_genome.bed -c > RE_sites_genome_N_windows_with_RE

#Count the number of RE>1 ((different cutoff instead can be used, like more than 2 reads, the top 10% and so on) in the window for the sample
awk '$4>1' RE_assgin_Sample.bedgraph | \
	intersectBed -a RE_sites_genome_N_windows -b stdin -c > RE_sites_genome_N_windows_with_Sample_RE 

#Use the files created in the R script.