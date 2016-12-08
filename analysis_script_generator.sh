#!/bin/bash
#Generate analysis script for every sample. The fastq files are generated from bam files in script "sort_and_bam2fastq.sh". Script content was written by my colleague Du Xiaohong.
for i in $(ls )
do
    echo "#!/bin/bash
cd /mnt/data2/Normal_people_mtDNA_data_from_limingkun/${i}
bwa aln -t 4 -f ${i:0:5}.R1.sai /home/database/mt_analysis/modified_genome/hg19_mt.fa /mnt/data2/Normal_people_mtDNA_data_from_limingkun/${i}/${i:0:5}_R1.fq
bwa aln -t 4 -f ${i:0:5}.R2.sai /home/database/mt_analysis/modified_genome/hg19_mt.fa /mnt/data2/Normal_people_mtDNA_data_from_limingkun/${i}/${i:0:5}_R2.fq
bwa sampe -r '@RG\tID:${i:0:5}\tLB:mtDNA\tPL:ILLUMINA\tSM:${i:0:5}' -f ${i:0:5}.raw.sam /home/database/mt_analysis/modified_genome/hg19_mt.fa ${i:0:5}.R1.sai ${i:0:5}.R2.sai /mnt/data2/Normal_people_mtDNA_data_from_limingkun/${i}/${i:0:5}_R1.fq /mnt/data2/Normal_people_mtDNA_data_from_limingkun/${i}/${i:0:5}_R2.fq
samtools view -bS ${i:0:5}.raw.sam -o ${i:0:5}.raw.bam 
rm ${i:0:5}.raw.sam 
java -jar /home/duxiaohong/software/picard/picard-tools-1.81/SortSam.jar I=${i:0:5}.raw.bam O=${i:0:5}.sort.picard.bam SO=coordinate VALIDATION_STRINGENCY=SILENT && echo "picard sort is finished " && date 
java -jar /home/duxiaohong/software/picard/picard-tools-1.81/BuildBamIndex.jar INPUT=${i:0:5}.sort.picard.bam OUTPUT=${i:0:5}.sort.picard.bai VALIDATION_STRINGENCY=SILENT && echo "picard index ${i:0:5}.sorted.bam is finished "&& date
java -jar /home/duxiaohong/software/picard/picard-tools-1.81/MarkDuplicates.jar REMOVE_DUPLICATES=true I=${i:0:5}.sort.picard.bam O=${i:0:5}.rmdup.bam M=${i:0:5}.dup.metric VALIDATION_STRINGENCY=SILENT && echo "picard mark duplicate is finished " && date 
java -jar /home/duxiaohong/software/picard/picard-tools-1.81/BuildBamIndex.jar INPUT=${i:0:5}.rmdup.bam  OUTPUT=${i:0:5}.rmdup.bai VALIDATION_STRINGENCY=SILENT && echo "picard index ${i:0:5}.sorted.bam is finished "&& date
#perl /home/duxiaohong/program/mt_analysis/double_barcode/diff_alignment_compare/fetch_mt.pl ${i:0:5}.rmdup.bam |samtools view -bS - -o ${i:0:5}.mt.bam 
#java -jar /home/duxiaohong/software/picard/picard-tools-1.81/BuildBamIndex.jar INPUT=${i:0:5}.rmdup.bam  OUTPUT=${i:0:5}.mt.bai VALIDATION_STRINGENCY=SILENT && echo "picard index ${i:0:5}.sorted.bam is finished "&& date
java -Xmx4g -jar /home/duxiaohong/software/GATK/GenomeAnalysisTK.jar -R /home/database/mt_analysis/modified_genome/hg19_mt.fa -T RealignerTargetCreator -o ${i:0:5}.realn.intervals -I ${i:0:5}.rmdup.bam && echo " GATK RealignerTargetCreator  is finished " && date 
java -Xmx4g -jar /home/duxiaohong/software/GATK/GenomeAnalysisTK.jar -R /home/database/mt_analysis/modified_genome/hg19_mt.fa -T IndelRealigner --maxReadsForRealignment 100000 -targetIntervals ${i:0:5}.realn.intervals -o ${i:0:5}.realn.bam -I ${i:0:5}.rmdup.bam &&  echo " GATK IndelRealigner is finished " && date 
#java -Xmx4g -jar /home/duxiaohong/software/GATK/GenomeAnalysisTK.jar -R /home/database/mt_analysis/modified_genome/human_mtDNA.fasta -T DepthOfCoverage --omitIntervalStatistics --omitLocusTable --minMappingQuality 0 -I ${i:0:5}.realn.bam -o ${i:0:5}.depth && echo " GATK depth is finished " && date 
samtools index ${i:0:5}.realn.bam 
perl /home/duxiaohong/program/mt_analysis/double_barcode/diff_alignment_compare/fetch_mt.pl ${i:0:5}.realn.bam |samtools view -bS - -o ${i:0:5}.mt.bam 
samtools index ${i:0:5}.mt.bam 
samtools depth -Q 0 ${i:0:5}.mt.bam > ${i:0:5}.depth
samtools view ${i:0:5}.mt.bam > ${i:0:5}.mt.sam 
python /home/duxiaohong/software/DREEP/sam_to_pileup.py ${i:0:5}.mt.sam /home/database/mt_analysis/modified_genome/human_mtDNA.fasta > ${i:0:5}.old.pileup 
perl /home/duxiaohong/software/DREEP/pileup_appended.pl ${i:0:5}.old.pileup > ${i:0:5}.pileup
perl /home/duxiaohong/software/DREEP/filter_and_summary_v2.pl -i ${i:0:5}.pileup -d 3 -q 30 -m 20 -r 118 -s 5  > ${i:0:5}.ssp
rm ${i:0:5}.old.pileup
rm ${i:0:5}.mt.sam
gzip -f ${i:0:5}.pileup
    " > ${i}/${i:0:5}_analysis.sh
done
