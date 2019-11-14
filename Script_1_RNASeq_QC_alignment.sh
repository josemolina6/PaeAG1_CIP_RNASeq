# 	Transcriptomic analysis of response of Pseudomonas aeruginosa AG1 after exposure to Ciprofloxacin

echo "SCRIPT FOR BIOINFORMATICS ANALYSIS"
echo "----------------------------------------------"
echo "Implemented by Jose Arturo Molina Mora"
echo "University of Costa Rica"
echo "----------------------------------------------"


echo "----------------------------------------------"
echo "SCRIPT 1: QUALILY CONTROL OF RAW DATA"
echo "----------------------------------------------"

# SAMPLES: A (TIME 0 h), B (TIME 2.5 h) AND C (5 h), by triplicates (numbers 1 to 3)


mkdir fastqc
fastqc -o fastqc *.fq.gz

mkdir fastqscreen
fastq_screen --aligner bwa *.fq.gz --outdir fastqscreen

multiqc ./fastqc/. -o fastqc
multiqc ./fastqscreen/. -o fastqscreen


# TRIMMING: ANALYSIS SHOWED FOR THE FIRST SAMPLE


for i in A1 A2 A3 B1 B2 B3 C1 C2 C3
do

trimmomatic PE ${i}_1P.fq.gz ${i}_2P.fq.gz -baseout ${i}_trim.fq.gz SLIDINGWINDOW:4:30 MINLEN:35 

echo "${i} FINISHED"
done

# FILTERING ADAPTS AND rRNA

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3
do

sh bbduk.sh in=${i}_trim_1P.fq.gz in2=${i}_trim_2P.fq.gz out=${i}_trim-fil_1P.fq.gz out2=${i}_trim-fil_2P.fq.gz outm=${i}_1P_OUT.fq.gz outm2=${i}_2P_OUT.fq.gz ref=adapters.fa,16S_rRNA_bact_archaea.fasta,23S_rRNA_bact_archaea.fasta stats=${i}_Stat.txt refstats=${i}_refstats.txt

echo "${i} FINISHED"
done

# FASTQC, FASTQ SCREEN AND MULTIQC WERE RUN AGAIN


echo "----------------------------------------------"
echo "SCRIPT 2: EDGE-PRO: mapping reads to PaeAG1 genome"
echo "----------------------------------------------"

# RUNNING EDGE-PRO FOR ALL SAMPLES, SHOWN FOR THE FIRST SAMPLE


for i in A1 A2 A3 B1 B2 B3 C1 C2 C3
do

perl EDGE_pro_v1.3.1/edge.pl -g PaeAG1.fasta \
        -p PaeAG1.ptt \
        -r PaeAG1.rnt \
        -u ${i}_trim-fil_1P.fq.gz \
        -o ${i}_Epro.out \
        -s EDGE_pro_v1.3.1/ \
        -v ${i}_trim-fil_2P.fq.gz \
        -t 8
echo "${i} FINISHED"
done

mv *.rpkm_0 ./fordeseq/raw_counts.txt

#Creating tables for DESEQ2:

perl EDGE_pro_v1.3.1/additionalScripts/edgeToDeseq.perl *.txt

# CONTINUE WITH SCRIPT 2 DESEQ2 (IN R SOFTWARE)

echo "----------------------------------------------"
echo "SCRIPT 3: MAPPING QUALITY CONTROL "
echo "----------------------------------------------" 

# BOWTIE

# USE ALIGNMENT FILES OF EDGE-PRO, OR DO ALIGNMENT USING BOWTIE2 TO GENERATE BAM FILES DIRECTLY TO EVALUATE A SIMILAR FILE ALIGMENT

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3
do
bowtie2 --local --no-unal -x PaeAG1.fasta -q -1 ${i}_trim-fil_1P.fq.gz,${i}_trim-fil_2P.fq.gz | samtools sort > ${i}_sort.bam
echo "${i} FINISHED"
done

# RUNNING QUALIMAP-RNASEQ

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3
do
qualimap rnaseq -bam ${i}_sort.bam -gtf AG1_ann_genome/PaeAG1.gtf -outdir qualimap_rnaseq/${i}_rnaseq -outformat PDF:HTML -pe --java-mem-size=10G
echo "${i} FINISHED"
done

multiqc ./qualimap_rnaseq/. -o qualimap_rnaseq


#RUNNING RSeQC

#gene_body_coverage, input is the folder with all bam files
geneBody_coverage.py -r AG1_ann_genome/PaeAG1.bed12 -i bam_files -o All_samples

#TIN (transcript integrity number)
tin.py -r AG1_ann_genome/PaeAG1.bed12 -i -i bam_files 


# END

