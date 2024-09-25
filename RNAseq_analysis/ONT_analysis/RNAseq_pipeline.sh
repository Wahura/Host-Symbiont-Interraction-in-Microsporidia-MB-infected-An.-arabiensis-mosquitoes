#!/bin/bash

##pycoQC
pycoQC -f sequencing_summary_FAR95474_39d7de8f.txt sequencing_summary_FAR95474_39d7de8f.txt -o pycoQC_output.html

##porechop

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/fastq_pass"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/output_reads"

for n in {01..12};
do mkdir -p "${output_dir}/barcode$n"

for i in ${input_dir}/barcode$n/*.fastq.gz;
do
base_name=$(basename "$i" .fastq.gz);

porechop -i "$i" --format fastq.gz -o "${output_dir}/barcode$n/${base_name}_porechopped.fastq.gz"

done
done

##concatenating files

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/output_reads"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/sortmerna_outputs"

mkdir -p "${output_dir}"

for i in "${input_dir}"/*/;
do
sub_dir_name=$(basename "$i")

output_file="${output_dir}/${sub_dir_name}_concatenated.fastq.gz"

cat "$i"/*.fastq.gz > "$output_file"

done

##sortmeRNA

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/sortmerna_outputs"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_input"
database_dir="/home/jwahura/lustre/RNAseq/data/72_PN/sortmerna_outputs/sortmerna_rRNA_databases"

mkdir -p ${output_dir}

for i in "${input_dir}"/*_concatenated.fastq.gz;
do
base_name=$(basename $i _concatenated.fastq.gz)

#inputs
input="${input_dir}/${base_name}_concatenated.fastq.gz"
db1="${database_dir}/rfam-5.8s-database-id98.fasta"
db2="${database_dir}/rfam-5s-database-id98.fasta"
db3="${database_dir}/silva-arc-16s-id95.fasta"
db4="${database_dir}/silva-arc-23s-id98.fasta"
db5="${database_dir}/silva-bac-16s-id90.fasta"
db6="${database_dir}/silva-bac-23s-id98.fasta"
db7="${database_dir}/silva-euk-18s-id95.fasta"
db8="${database_dir}/silva-euk-28s-id98.fasta"

#outputs
out1="${output_dir}/${base_name}.aligned.porechopped"
out2="${output_dir}/${base_name}.unaligned.porechopped"
work="${input_dir}"

sortmerna --ref $db1 --ref $db2 --ref $db3 --ref $db4 --ref $db5 --ref $db6 --ref $db7 --ref $db8 --reads $input --aligned $out1 --fastx --other $out2 --workdir $work

rm -r "${work}/kvdb/*"

done

##unzipping the .gz files

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_input"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_input"

for i in "${input_dir}"/*.unaligned.porechopped.fq.gz;
do

gunzip "$i"

done

##isonclust

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_input"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_output"

for i in "${input_dir}"/*.unaligned.porechopped.fq; 
do
 
base_name=$(basename "$i" .unaligned.porechopped.fq)

mkdir -p "${output_dir}/${base_name}"

isONclust --ont --fastq "$i" --outfolder "${output_dir}/${base_name}"

done

##isonclust writing

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_output"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isoncorrect_input"
files="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_input"

for i in "${files}"/*.unaligned.porechopped.fq;

do

base_name=$(basename "$i" .unaligned.porechopped.fq)

mkdir -p "${output_dir}/${base_name}"

isONclust write_fastq --clusters "${input_dir}/${base_name}/final_clusters.tsv" --fastq "${files}/${base_name}.unaligned.porechopped.fq" --outfolder "${output_dir}/${base_name}" --N 1

done

##isoncorrect

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isoncorrect_input"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isoncorrect_output"
files="/home/jwahura/lustre/RNAseq/data/72_PN/isonclust_input"

for i in "${files}"/*.unaligned.porechopped.fq;

do

base_name=$(basename "$i" .unaligned.porechopped.fq)

mkdir -p "${output_dir}/${base_name}_correction"

run_isoncorrect --fastq "${input_dir}/${base_name}" --outfolder "${output_dir}/${base_name}" --exact_instance_limit 50 --max_seqs 1000 --k 9 --w 10 --xmin 14 --xmax 80 --T 0.1

done

##merging all corrected files into one

input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/isoncorrect_output"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_input"

mkdir -p "${output_dir}"

barcodes=("barcode01" "barcode02" "barcode03" "barcode04" "barcode05" "barcode06" "barcode07" "barcode08" "barcode09" "barcode10" "barcode11" "barcode12")

for barcode in ${barcodes[@]};
do

echo ${input_dir}/${barcode}/*/corrected_reads.fastq

cat ${input_dir}/${barcode}/*/corrected_reads.fastq >> ${output_dir}/${barcode}_all_corrected_reads.fq

done

##mapping using minimap
input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_input"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_output"

mkdir -p "${output_dir}"

for i in ${input_dir}/*_all_corrected_reads.fq

do 

base_name=$(basename "$i" _all_corrected_reads.fq)

ref="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_input/VectorBase-63_AarabiensisDONGOLA2021_AnnotatedTranscripts.fasta"


minimap2 -ax map-ont -k13 ${ref} "$i" > ${output_dir}/${base_name}_aln.sam 

done

##converting sam to bam files
input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_output"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_output"


for i in ${input_dir}/*_aln.sam; 
do
base_name=$(basename "$i" _aln.sam)
  
samtools view -S -b "$i" | samtools sort - | samtools view -h -F 4 - > "${output_dir}/${base_name}_bwa_mapped.bam"

done

##counting using nanocount
input_dir="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_output"
output_dir="/home/jwahura/lustre/RNAseq/data/72_PN/minimap_output"


for i in ${input_dir}/*_bwa_mapped.bam; 
do
base_name=$(basename "$i" _bwa_mapped.bam)
  
NanoCount -i "$i" -o "$output_dir/${base_name}_counts.tsv"

done



















