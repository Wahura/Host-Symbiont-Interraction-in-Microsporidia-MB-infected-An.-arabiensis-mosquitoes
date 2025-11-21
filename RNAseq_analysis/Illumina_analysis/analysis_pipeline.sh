#!/bin/bash

##fastqc
input_dir="/home/jwahura/lustre/RNAseq/data/aws_data/fastq_combined"
output_dir="/home/jwahura/lustre/RNAseq/data/aws_data/fastqc_output"

mkdir -p "${output_dir}"

fastqc -t 24 -o "${output_dir}" ${input_dir}/*.fq.gz

##multiqc
input_dir="/home/jwahura/lustre/RNAseq/data/aws_data/fastqc_output"
output_dir="/home/jwahura/lustre/RNAseq/data/aws_data/multiqc_report"

mkdir -p "${output_dir}"

multiqc "${input_dir}" -o "${output_dir}"

##sortmeRNA
input_dir="/home/jwahura/lustre/RNAseq/data/aws_data/fastq_combined"
output_dir="/home/jwahura/lustre/RNAseq/data/aws_data/hisat_input"
database_dir="/home/jwahura/lustre/RNAseq/data/aws_data/sortmerna_rRNA_databases"

mkdir -p ${output_dir}

for i in "${input_dir}"/*.fq.gz;
do
base_name=$(basename $i .fq.gz)

##inputs
input="${input_dir}/${base_name}.fq.gz"
db1="${database_dir}/rfam-5.8s-database-id98.fasta"
db2="${database_dir}/rfam-5s-database-id98.fasta"
db3="${database_dir}/silva-arc-16s-id95.fasta"
db4="${database_dir}/silva-arc-23s-id98.fasta"
db5="${database_dir}/silva-bac-16s-id90.fasta"
db6="${database_dir}/silva-bac-23s-id98.fasta"
db7="${database_dir}/silva-euk-18s-id95.fasta"
db8="${database_dir}/silva-euk-28s-id98.fasta"

##outputs
out1="${output_dir}/${base_name}.aligned"
out2="${output_dir}/${base_name}.unaligned"
work="${input_dir}"

sortmerna --ref $db1 --ref $db2 --ref $db3 --ref $db4 --ref $db5 --ref $db6 --ref $db7 --ref $db8 --reads $input --aligned $out1 --fastx --other $out2 --workdir $work

rm /home/jwahura/lustre/RNAseq/data/aws_data/fastq_combined/kvdb/*

done

##splitting the reads
input_dir="/home/jwahura/lustre/RNAseq/data/aws_data/hisat_input"
output_dir="/home/jwahura/lustre/RNAseq/data/aws_data/split_reads"
mkdir -p "${output_dir}"

for fq in "${input_dir}"/*.unaligned.fq.gz; do
    sample=$(basename "$fq" .unaligned.fq.gz)
    echo "Splitting $sample..."

    gunzip -c "$fq" | awk -v s="${sample}" -v out="${output_dir}" '
    BEGIN {
        r1 = out "/" s "_1.fq"
        r2 = out "/" s "_2.fq"
    }
    {
        header = $0
        getline seq
        getline plus
        getline qual

        if (header ~ /\/1$/) {
            print header >> r1
            print seq   >> r1
            print plus  >> r1
            print qual  >> r1
            fflush(r1)
        } else if (header ~ /\/2$/) {
            print header >> r2
            print seq   >> r2
            print plus  >> r2
            print qual  >> r2
            fflush(r2)
        }
    }'

    echo "$sample done."
done



##hisat2
input_dir="/home/jwahura/lustre/RNAseq/data/aws_data/split_reads"
output_dir="/home/jwahura/lustre/RNAseq/data/aws_data/hisat_mapped"

mkdir -p "${output_dir}"

for sample in $(ls ${input_dir}/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | xargs -n 1 basename); do
    hisat2 -p 24 \
        -x anopheles_index \
        -1 "${input_dir}/${sample}_R1.fastq.gz" \
        -2 "${input_dir}/${sample}_R2.fastq.gz" \
        -S "${output_dir}/${sample}.sam"
done

##samtools to convert sam files to bam files
input_dir="/home/jwahura/lustre/RNAseq/data/aws_data/hisat_mapped"
output_dir="/home/jwahura/lustre/RNAseq/data/aws_data/hisat_mapped"


for i in ${input_dir}/*.sam; 
do
base_name=$(basename "$i" .sam)
  
samtools view -S -b "$i" | samtools sort - | samtools view -h -F 4 - > "${output_dir}/${base_name}_hisat_mapped.bam"

done

##samtools to quantify the transcripts
input_dir="/home/jwahura/lustre/RNAseq/data/aws_data/transcriptome/hisat_mapped"
output_dir="/home/jwahura/lustre/RNAseq/data/aws_data/transcriptome/hisat_mapped"


for i in ${input_dir}/*.bam; 
do
base_name=$(basename "$i" .bam)
  
 samtools view -F 4 "$i" | cut -f 3 | sort | uniq -c > "${output_dir}/${base_name}_counts.tsv"
done






