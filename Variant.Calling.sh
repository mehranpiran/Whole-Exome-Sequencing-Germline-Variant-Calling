
############ Download fastq files using sratoolkit programme #####

Download the programme
$ wget (sratoolkit software download link)

uncompress the file.tar.gz
$ tar -zxvf sratoolkit.2.9.2_ubuntu_64.tar.gz


Run the programme for paired-end
chipseq data
$ ~/sratoolkit.2.9.2_ubuntu_64/bin/fastq-dump -I --split-files SRRXXXXXXX > treatment.fastq
input data
$ ~/sratoolkit.2.9.2_ubuntu_64/bin/fastq-dump -I --split-files SRRXXXXXXX > input.fastq


Run the programme for single-end
chipseq data
$ ~/sratoolkit.2.9.2_ubuntu_64/bin/fastq-dump -z SRRXXXXXXX > treatment.fastq
input data
$ ~/sratoolkit.2.9.2_ubuntu_64/bin/fastq-dump -z SRRXXXXXXX > input.fastq 


############ separate forward and revers reads #####

Use sratoolkit
OR
awk '0 == ((NR+4) % 8)*((NR+5) % 8)*((NR+6) % 8)*((NR+7) %8)' SRRXXXXX_sra.fastq > SRRXXXXX_sra_1.fastq  
awk '0 == (NR % 8)*((NR+1) % 8)*((NR+2) % 8)*((NR+3) %8)' SRRXXX_sra.fastq > SRRXXX_sra_2.fastq 


########### Quality control of the reads #####

$ wget (fastqc software download link)
$ gunzip fastqc_v0.11.7.zip
$ cd FastQC
$ chmod u+x fastqc
$ ~/apps/FastQC/fastqc ~/Data/*.fastq




########### Trim reads #####

$ wget (Trimmomatic software download link)
$ gunzip Trimmomatic-0.36.zip

For paired-end
$ java -jar ~/bin/trimmomatic-0.36.jar PE treatment_1.fastq treatment_2.fastq treatment_1_paired_trimmed.fq.gz treatment_1_unpaired.fq.gz treatment_2_paired_trimmed.fq.gz treatment_2_unpaired.fq.gz LEADING:10 TRAILING:10 MINLEN:25 HEADCROP:5 TAILCROP:5 SLIDINGWINDOW:4:20
or
$ for f in *_1.fastq; do java -jar ~/apps/Trimmomatic-0.36/trimmomatic-0.36.jar PE $f ${f/_1/_2} -baseout ${f/_treatment*/.fastq} HEADCROP:5 LEADING:10 TRAILING:10 MINLEN:25; done

For Single-end
$ java -jar ~/apps/Trimmomatic-0.36/trimmomatic-0.36.jar SE treatment.fastq treatment_paired_trimmed.fq.gz treatment_unpaired.fq.gz LEADING:10 TRAILING:10 MINLEN:25 HEADCROP:5 TAILCROP:5 SLIDINGWINDOW:4:20
or
for f in *.fastq; do java -jar ~/apps/Trimmomatic-0.36/trimmomatic-0.36.jar SE $f HEADCROP:5 LEADING:10 TRAILING:10 MINLEN:25; done



############################                 #################################
############################  Preprocessing  #################################
############################                 #################################


1) Alignment with BWA

#change to reference sequence file directory
cd ~/Science/BWA.index.alignment

1.1) Make the reference bwa Index
bwa index -a bwtsw hg38.fa
#this function creates multiple files with different extensions. It is a highly energy-consuming process.
#configured files are amb, ann, bwt, pac, sa.

1.2) start the alignment
bwa mem -M -R '@RG\tID:group1\tSM:skbr3\tPL:illumina' hg38.fa SRR8597579_1_paired_trimmed.fq.gz SRR8597579_2_paired_trimmed.fq.gz > SRR8597579.sam
*************************
-M for split alingment. Split reads: When one portion of an NGS read map to one location and other portion of the same read map to a different location of a genome. When the read both the portions is of equal length, this is called a balanced split.

-R '@RG\tID:group1\tSM:skbr3\tPL:illumina' for header in sam file. -R <readgroup_info> Provide the readgroup as a string. The read group information is key for downstream GATK functionality. The GATK will not work without a read group tag.
*************************

2) Sort sam and convert to bam
java -jar ~/apps/picard.jar SortSam I=sra_data.sam O=sra_data.bam SORT_ORDER=coordinate

#java -jar build/libs/picard.jar SortSam -I testdata/picard/sam/namesorted.test.sam -SORT_ORDER coordinate -O sorted.sam


# To view the bam files "samtools" is used
samtools view sra_data.bam | head -5


3) Alignment Metrics
#Let’s compute some statistics to see how well our reads aligned to the reference genome. We’ll use samtools flagstat for this.
samtools flagstat sra_data.sam


4) Mark Duplicates
# We use Picard Tools to locate and tag duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Note that this step does not remove the duplicate reads, but rather flags them as such in the read’s SAM record.

java -jar ~/apps/picard.jar MarkDuplicates INPUT=sra_data.bam OUTPUT=sra_data_dup.bam METRICS_FILE=dupmetrics.txt

*************************
link for explains about Picard
http://broadinstitute.github.io/picard/command-line-overview.html
*************************

5) Prepare reference dictionary, fasta index, and bam index 

5.1) Make the fasta file index
samtools faidx hg38.fa
#It creates hg38.fa.fai file


5.2) make reference dictionary
java -jar ~/apps/picard.jar CreateSequenceDictionary REFERENCE=~/Science/BWA.index.alignment/hg38.fa OUTPUT=hg38.dict
# It generates hg38.dic file


5.3) Index Bam file:Generates a BAM index ".bai" file. This tool creates an index file for the input BAM that allows fast look-up of data in a BAM file, lke an index on a database. Note that this tool cannot be run on SAM files, and that the input BAM file must be sorted in coordinate order.

samtools index sra_data_dup.bam
or
java -jar ~/apps/picard.jar BuildBamIndex I=sra_data.bam o=sra_data.bam.bai

# It generates sra_data_dup.bam.bai file.



6) Base recalibration

http://data.broadinstitute.org/igvdata/annotations/hg38/snp142.snp.gz

java -jar $GATK -T BaseRecalibrator -R reference/chr19.fa -I chr19_dedup.bam -knownSites latest_dbsnp.vcf -o recal_data.table

6.1) java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator -R hg38.fa -I sra_data_dup.bam -knownSites 1000GENOMES-phase_3.vcf.gz -o recal_data.table

6.2) java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T PrintReads -R hg38.fa -I sra_data_dup.bam -BQSR recal_data.table -o recal_reads.bam



################## Variant Discovery ##################

1) Calling variants

1.1) java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R hg38.fa -I recal_reads.bam -o raw_variants.vcf

1.2) separate SNPs from InDels

java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SelectVariants -R hg38.fa -V raw_variants.vcf -selectType SNP -o raw_snps.vcf

java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SelectVariants -R hg38.fa -V raw_variants.vcf -selectType INDEL -o raw_indels.vcf



2) Filter Variants


2.1) SNPs

java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R hg38.fa -V raw_snps.vcf -filterName "QD_filter" -filter "QD<2.0" -filterName "FS_filter" -filter "FS>60.0" -filterName "MQ_filter" -filter "MQ<40.0" -filterName "SOR_filter" -filter "SOR>10.0" -o filtered_snps.vcf

OR

java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration \
-R hg38.fa \
-V raw_snps.vcf \
-filterName "QD_filter" \
-filter "QD<2.0" \
-filterName "FS_filter" \
-filter "FS>60.0" \
-filterName "MQ_filter" \
-filter "MQ<40.0" \
-filterName "SOR_filter" \
-filter "SOR>10.0" \
-o filtered_snps.vcf


2.2) InDels

java -jar ~/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R hg38.fa -V raw_indels.vcf -filterName "QD_filter" -filter "QD<2.0" -filterName "FS_filter" -filter "FS>200.0" -filterName "SOR_filter" -filter "SOR>10.0" -o filtered_indels.vcf

OR

java -jar $GATK_JAR -T VariantFiltration \
-R GCF_000001405.33_GRCh38.p7_chr20_genomic.fna \
-V raw_indels.vcf \
-filterName "QD_filter" \
-filter "QD<2.0" \
-filterName "FS_filter" \
-filter "FS>200.0" \
-filterName "SOR_filter" \
-filter "SOR>10.0" \
-o filtered_indels.vcf




3) Annotation

java -jar snpEff.jar -v snpeff_db filtered_snps_final.vcf > filtered_snps_final.ann.vcf

java -jar ~/apps/snpEff/snpEff.jar -v ~/Science/Variant Discovery/00-All.vcf_3.gz -ann filtered_snps.vcf > filtered_snps_ann.vcf

