# 16Samplicons

Following Moving pictures tutorial for Qiime2 and IMR's workflow https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-Helper-2-Marker-gene-workflow, but with modifications for paired-end reads.

All analyses are completed on OSC unless stated otherwise: see PGG-OSC-Intro for instructions on getting started with OSC and transferring files in.

My files are in /users/PAS3057/qfaber/quelccaya/seqs/ and everything is being run in /users/PAS3057/qfaber/qiime2_test

Note to self: eventually add in our own practice dataset with the primers we use! How can i have blocks of code within the text like the IMR tutorial has?? and also like not in a readme I guess. Apply to other module too!
Eventually add in things to help decipher what a good quality score is, etc. The interactive visuals are great so use them!!!


##Installing Qiime2

<pre> module load miniconda3/24.1.2-py310 </pre>

<pre> conda env create \
  --prefix /users/PAS3057/qfaberconda3/envs/amplicon-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml </pre>

This will take longer to install than normal programs with conda due to the amount of packages required.

<pre> conda activate /users/PAS3057/qfaber/miniconda3/envs/qiime2 </pre>

Test Install 
  
<pre> qiime info </pre>

May need to install additional programs

I was missing a program in R so I did the following:

R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData", ask = FALSE)

library(GenomeInfoDbData)
library(phyloseq)

q()

qiime info

## Basic quality check of sequences

##Install fastqc and multiqc

conda create --prefix /users/PAS3057/qfaber/miniconda3/envs/qualitycheck -c bioconda fastqc multiqc
conda activate /users/PAS3057/qfaber/miniconda3/envs/qualitycheck
mkdir fastqc_out
fastqc -t 1 /users/PAS3057/qfaber/quelccaya/seqs/*.fastq.gz -o fastqc_out
multiqc fastqc_out --filename multiqc_report.html

Look at overall quality of your reads.

##Prepare metadata

Data should be in a .tsv file, which can be created in excel. First column should be sample_id and all ids should be unique.

## importing data
This assumes that file names are in the default format provided by AMR that have not been renamed. For more details on importing other files types see https://amplicon-docs.qiime2.org/en/stable/how-to-guides/how-to-import.html 

mkdir reads_qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /users/PAS3057/qfaber/quelccaya/seqs/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path reads_qza/reads.qza


qiime demux summarize \
   --i-data reads_qza/reads.qza \
   --o-visualization reads_qza/reads_summary.qzv

Upload file to https://view.qiime2.org/ to look at interactive quality plots

## Trim primers with cutadapt

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences reads_qza/reads.qza \
  --p-cores 8 \
  --p-anywhere-f GTGYCAGCMGCCGCGGTAA \
  --p-anywhere-r GGACTACNVGGGTWTCTAAT \
  --p-error-rate 0.1 \
  --p-match-read-wildcards true \
  --p-match-adapter-wildcards true \
  --o-trimmed-sequences reads_qza/reads_trimmed.qza 


Note: this is assuming 16S V4 primers, but will change if other primers are used!
Original code trimmed absolutely everything

Take a look at the sequence quality once again:

qiime demux summarize \
   --i-data reads_qza/reads_trimmed.qza \
   --o-visualization reads_qza/reads_trimmed.qzv

qiime tools export \
  --input-path reads_qza/reads_trimmed.qza \
  --output-path reads_fastq

mkdir -p reads_fastq/merged reads_fastq/unmerged



##Denoising

qiime tools export \
  --input-path reads_qza/reads_trimmed.qza \
  --output-path reads_fastq

conda create --prefix /users/PAS3057/qfaber/miniconda3/envs/bbtools -c bioconda fastp
conda activate /users/PAS3057/qfaber/miniconda3/envs/fastp

mkdir -p reads_fastq/merged reads_fastq/unmerged

for fwd in reads_fastq/*_R1_001.fastq.gz
do
    sample=$(basename "$fwd" _R1_001.fastq.gz)
    rev="reads_fastq/${sample}_R2_001.fastq.gz"

    fastp \
      -i "$fwd" \
      -I "$rev" \
      -m \
      --merged_out "reads_fastq/merged/${sample}_merged.fastq.gz" \
      --unpaired1 "reads_fastq/unmerged/${sample}_unmerged_R1.fastq.gz" \
      --unpaired2 "reads_fastq/unmerged/${sample}_unmerged_R2.fastq.gz" \
      --length_required 200 \
      --overlap_diff_percent_limit 10 \
      --overlap_len_require 30 \
      --thread 8
done
Note: the length required will be different depending on whether you are doing V4 or V4V5.

conda activate /users/PAS3057/qfaber/miniconda3/envs/qiime2

cd reads_fastq/merged

# Create header for manifest
echo "sample-id,absolute-filepath,direction" > ~/qiime2_test/reads_qza/manifest.csv

# Loop through all merged FASTQs
for f in *_merged.fastq.gz; do
    sample=$(echo $f | sed 's/_S[0-9]\+_L001_merged.fastq.gz//')  # remove suffix to get sample ID
    abspath=$(realpath "$f")  # get full path
    echo "${sample},${abspath},forward" >> ~/qiime2_test/reads_qza/manifest.csv
done

cat reads_qza/manifest.csv | tr ',' '\t' > reads_qza/manifest.tsv

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path reads_qza/manifest.tsv \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path reads_qza/reads_merged.qza

qiime quality-filter q-score \
   --i-demux reads_qza/reads_merged.qza \
   --o-filter-stats filt_stats.qza \
   --o-filtered-sequences reads_qza/reads_mered_filt.qza

qiime demux summarize \
   --i-data reads_qza/reads_mered_filt.qza \
   --o-visualization reads_qza/reads_mered_filt_summary.qzv
Adjust for number of cores you are using and length you want to trim to depending on primer set

qiime deblur denoise-16S \
  --i-demultiplexed-seqs reads_qza/reads_mered_filt.qza \
  --p-trim-length 253 \
  --p-sample-stats \
  --p-jobs-to-start 8 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --o-stats deblur-stats.qza

Look at deblur statistics:

qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv

mv rep-seqs-deblur.qza rep-seqs.qza
mv table-deblur.qza table.qza


Get sequences of all:

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv


##taxonomic classification

wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza

Note: run on 8 cores since it takes a while and request 100G of RAM
qiime feature-classifier classify-sklearn \
   --i-reads rep-seqs.qza \
   --i-classifier silva-138-99-nb-classifier.qza \
   --p-read-orientation same \
   --p-n-jobs 8 \
   --output-dir taxa


EXPORT ALL TO R:

qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy

This will save tsv file under exported-taxonomy

Enter rep-seqs.qzv to view.qiime2.org and download fasta file

qiime tools export \
  --input-path deblur_output/table.qza \
  --output-path exported-feature-table
biom convert \
  -i exported-feature-table/feature-table.biom \
  -o feature-table.tsv \
  --to-tsv

This will give tsv with number of counts per sample
