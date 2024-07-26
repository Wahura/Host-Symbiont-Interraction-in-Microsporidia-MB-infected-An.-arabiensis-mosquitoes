#!/bin/bash

#Importing data into qiime2
qiime tools import \ 
--type 'SampleData[PairedEndSequencesWithQuality]' \ 
--input-path data \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path demux_paired_end.qza

# Visulize reads (QA)
qiime demux summarize \
--i-data demux_paired_end.qza \
--o-visualization demux_paired_end.qzv

#denoising
qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux_paired_end.qza \
--p-trunc-len-f 260 \
--p-trunc-len-r 200 \
--p-trim-left-f 20 \
--p-trim-left-r 8 \
--o-table feature_table.qza \
--o-representative-sequences rep_seqs.qza \
--o-denoising-stats denoising_stats.qza

#visualizing the feature table
qiime feature-table summarize \
--i-table feature_table.qza \
--o-visualization feature_table.qzv

#visualizing the denoising stats 
qiime metadata tabulate \
--m-input-file denoising_stats.qza \
--o-visualization denoising_stats.qzv

#visualizing the rep_seqs
qiime feature-table tabulate-seqs \
--i-data rep_seqs.qza \
--o-visualization rep_seqs.qzv

# Taxonomic classification
qiime feature-classifier classify-sklearn \
--i-classifier silva138_AB_V3-V4_classifier.qza \
--p-n-jobs 1 \
--i-reads rep_seqs.qza \
--o-classification taxonomy.qza

#visualizing the taxonomy file
qiime metadata tabulate \ 
--m-input-file taxonomy.qza \ 
--o-visualization taxonomy.qzv 

#Exporting the feature table
qiime tools export \
  --input-path feature_table.qza \
  --output-path exported-feature_table

#converting feature table biom file to tsv
biom convert -i exported-feature-table/feature-table.biom -o 
feature-table.tsv --to-tsv

#converting rep-seq fasta file to tsv
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print "#"$1,$2}' sequences.fasta > 
mosquito_sequence.tsv

#Removing the # in the sequence tsv file
sed 's/#//' mosquito_sequence.tsv >modified_mosquito_sequence.tsv

#Running blast
blastn -db 16SMicrobial -query mosquito_16S_sequences.fasta 
-max_target_seqs 1 -out mosquito_16S.csv -outfmt '6 qseqid sseqid evalue 
bitscore sgi sacc staxids sskingdoms sscinames stitle'

#Running piecrust pipeline
#run piecrust pipeline
picrust2_pipeline.py -s picrust2_seq.fasta -i picrust2_feature_table.biom 
-o picrust2_out -p 20
