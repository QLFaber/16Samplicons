# 16Samplicons

Following Moving pictures tutorial for Qiime2 and IMR's workflow https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-Helper-2-Marker-gene-workflow, but with modifications for paired-end reads.

All analyses are completed on OSC unless stated otherwise: see PGG-OSC-Intro for instructions on getting started with OSC and transferring files in.

My files are in /users/PAS3057/qfaber/quelccaya/seqs/ and everything is being run in /users/PAS3057/qfaber/try2 CHANGE PATH

Note to self: eventually add in our own practice dataset with the primers we use! How can i have blocks of code within the text like the IMR tutorial has?? and also like not in a readme I guess. Apply to other module too!
Eventually add in things to help decipher what a good quality score is, etc. The interactive visuals are great so use them!!!


## Installing Qiime2

<pre> module load miniconda3/24.1.2-py310 </pre>

<pre> conda env create \
  --prefix /users/PAS3057/qfaberconda3/envs/amplicon-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml </pre>

This will take longer to install than normal programs with conda due to the amount of packages required.

<pre> conda activate /users/PAS3057/qfaber/miniconda3/envs/qiime2 </pre>

Test Install 
  
<pre> qiime info </pre>

May need to install additional programs: I was missing a program in R so I did the following:

<pre> R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData", ask = FALSE)

library(GenomeInfoDbData)
library(phyloseq)

q()
</pre>

Now Test Install again
<pre> qiime info </pre>

## Basic quality check of sequences

Install fastqc and multiqc

<pre>conda create --prefix /users/PAS3057/qfaber/miniconda3/envs/qualitycheck -c bioconda fastqc multiqc
conda activate /users/PAS3057/qfaber/miniconda3/envs/qualitycheck
mkdir fastqc_out
fastqc -t 1 /users/PAS3057/qfaber/quelccaya/seqs/*.fastq.gz -o fastqc_out
multiqc fastqc_out --filename multiqc_report.html </pre>

Look at overall quality of your reads.

## Importing data
This assumes that file names are in the default format provided by AMR that have not been renamed. For more details on importing other files types see https://amplicon-docs.qiime2.org/en/stable/how-to-guides/how-to-import.html 

<pre> mkdir reads_qza 
conda activate /users/PAS3057/qfaber/miniconda3/envs/qiime2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /users/PAS3057/qfaber/quelccaya/seqs/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path reads_qza/reads.qza

qiime demux summarize \
   --i-data reads_qza/reads.qza \
   --o-visualization reads_qza/reads_summary.qzv

  </pre>

Upload file to https://view.qiime2.org/ to look at interactive quality plots

## Trim primers with cutadapt

<pre>
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences reads_qza/reads.qza \
  --p-cores 8 \
  --p-anywhere-f GTGYCAGCMGCCGCGGTAA \
  --p-anywhere-r GGACTACNVGGGTWTCTAAT \
  --p-error-rate 0.1 \
  --p-match-read-wildcards true \
  --p-match-adapter-wildcards true \
  --o-trimmed-sequences reads_qza/reads_trimmed.qza 
</pre>

Note: this is assuming 16S V4 primers, but will change slightly if other primers are used!
Original code trimmed absolutely everything

Take a look at the sequence quality once again:

<pre>
qiime demux summarize \
   --i-data reads_qza/reads_trimmed.qza \
   --o-visualization reads_qza/reads_trimmed.qzv

qiime tools export \
  --input-path reads_qza/reads_trimmed.qza \
  --output-path reads_fastq
</pre>



## Merging in Fastp

<pre>
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
</pre>
Note: the length required will be different depending on whether you are doing V4 or V4V5.


## Import back into QIIME2
<pre>
conda activate /users/PAS3057/qfaber/miniconda3/envs/qiime2
cd reads_fastq/merged
echo "sample-id,absolute-filepath,direction" > /users/PAS3057/qfaber/try2/reads_qza/manifest.csv

for f in *_merged.fastq.gz; do
    sample=$(echo $f | sed 's/_S[0-9]\+_L001_merged.fastq.gz//')  # remove suffix to get sample ID
    abspath=$(realpath "$f")  # get full path
    echo "${sample},${abspath},forward" >> /users/PAS3057/qfaber/try2/reads_qza/manifest.csv
done

cd /users/PAS3057/qfaber/try2
  
cat reads_qza/manifest.csv | tr ',' '\t' > reads_qza/manifest.tsv

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path reads_qza/manifest.tsv \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path reads_qza/reads_merged.qza
</pre>

## Denoising with deblur
<pre>
qiime quality-filter q-score \
   --i-demux reads_qza/reads_merged.qza \
   --o-filter-stats filt_stats.qza \
   --o-filtered-sequences reads_qza/reads_merged_filt.qza

qiime demux summarize \
   --i-data reads_qza/reads_merged_filt.qza \
   --o-visualization reads_qza/reads_merged_filt_summary.qzv
</pre>


Adjust for number of cores you are using and length you want to trim to depending on primer set. I recommend doing this in a job with multiple cores:

<pre>
qiime deblur denoise-16S \
  --i-demultiplexed-seqs reads_qza/reads_merged_filt.qza \
  --p-trim-length 253 \
  --p-sample-stats \
  --p-jobs-to-start 8 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --o-stats deblur-stats.qza
</pre>


Look at deblur statistics:
<pre>
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv

mv rep-seqs-deblur.qza rep-seqs.qza
mv table-deblur.qza table.qza

</pre>




## Taxonomic classification

<pre>
wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
</pre>


Note: run on 8 cores since it takes a while and request 100G of RAM:
<pre>
qiime feature-classifier classify-sklearn \
   --i-reads rep-seqs.qza \
   --i-classifier silva-138-99-nb-classifier.qza \
   --p-read-orientation same \
   --p-n-jobs 8 \
   --output-dir taxa
</pre>

## Export for analysis in R

Input this into https://view.qiime2.org/ then download fasta file of all seqs
<pre>
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
</pre>

This will save tsv file under exported-taxonomy:
<pre>
qiime tools export \
  --input-path taxa/classification.qza \
  --output-path exported-taxonomy
</pre>

This will give tsv with number of counts per sample:
<pre>
qiime tools export \
  --input-path deblur_output/table.qza \
  --output-path exported-feature-table
biom convert \
  -i exported-feature-table/feature-table.biom \
  -o feature-table.tsv \
  --to-tsv
</pre>
