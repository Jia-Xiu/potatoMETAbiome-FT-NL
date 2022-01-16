#!/bin/sh
#SBATCH --time=50:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --job-name=import_seq
#SBATCH -o import_%j.out
#SBATCH -e import_%j.error
#SBATCH --mem 16G
#SBATCH --partition=gelifes
#SBATCH --mail-user=x.jia@rug.nl
#SBATCH --mail-type=ALL

# this script created by Xiu Jia (xibeihenai@gmail.com) on 06-01-2021
# updated by Xiu Jia on 16-01-2022
# it follows the scripts used in potatoMETAbiome project made by Stefanie Vink 
# sequences are analyzed using QIIME2 on Peregrine HPC (the University of Groingen). To access Peregrine, you can used MobaXterm in Windows and Terminal on MacOS
# in QIIME2 the qza files are the data files and the qzv are the visualization files
# all the visualization file in *.qzv format can be see at https://view.qiime2.org/

# load QIIME2 module. For more details about this version QIIME2, go https://docs.qiime2.org/2020.8/
module load QIIME2/2020.8


# define the name of the dataset
SOURCE="PMB-FT-NL" 
DB="/data/p278113/qiime/QIIME2/silva-138-99-515-806-nb-classifier.qza"

### first set the right working directory dataset and upload metadata and manifest file in correct formats
$PWD



### Using google sheet add-on, keemei, to check if the metadata file is correct


## import the sequence files # this takes a few mins depending on file size
## input types: forward, reverse and barcode fastq files
## output type: PairedEndSequencesWithQuality

  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest-$SOURCE.csv \
    --output-path demux-PE-$SOURCE.qza\
    --input-format PairedEndFastqManifestPhred33
    
# visualize imported sequences (eg. quality and length)
  qiime demux summarize \
    --i-data demux-PE-$SOURCE.qza \
    --o-visualization demux-PE-$SOURCE.qzv

   
# remove primer sequences using cutadapt - this will also remove any adapter sequences 
# if there are empty sequences in a file, you need to use cutadapt outside QIIME with -m 1 flag (removes any sequences
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux-PE-$SOURCE.qza \
    --p-cores 16 \
    --p-front-f GTGCCAGCMGCCGCGGTAA \
    --p-front-r GGACTACHVGGGTWTCTAAT \
    --o-trimmed-sequences demux-PE-trim-$SOURCE.qza


# visualize trimmed sequences
 qiime demux summarize \
    --i-data demux-PE-trim-$SOURCE.qza \
    --o-visualization demux-PE-trim-$SOURCE.qzv
  

## denoise and derep using dada2
## set time to 24 hours --ntasks=16, mem=50 GB <- including chimera checking ## was done in about 3 1/2 hours
## for chimera-method none you must use QIIME2 2018.2 or higher - options: consensus, pooled or none
## trim and trunc values should be determined from demuxPE.qzv graphs
 
## input type: SampleData[PairedEndSequencesWithQuality]
## output types: FeatureTable[Frequency] and FeatureData[Sequence]
 
 
 
# gh-cutadapt-trimmed input file: trim left-f and left-r: 0, trunc-f 231, trunc-r 230 (0.4% of seqs were shorter than this length!)

## denoise and derep using dada2
## set time to 24 hours --ntasks=16, mem=50 GB <

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux-PE-trim-$SOURCE.qza \
    --o-table table-$SOURCE.qza \
    --o-representative-sequences rep-seqs-$SOURCE.qza \
    --o-denoising-stats stats-dada2-$SOURCE.qza \
    --p-chimera-method pooled \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-trunc-len-f 230 \
    --p-trunc-len-r 220 \
    --p-n-threads 6 \
    --verbose


# FeatureTable and FeatureData summarize
qiime feature-table summarize \
  --i-table table-$SOURCE.qza \
  --o-visualization table-$SOURCE.qzv \
  --m-sample-metadata-file metadata-$SOURCE.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-$SOURCE.qza \
  --o-visualization rep-seqs-$SOURCE.qzv

qiime metadata tabulate \
  --m-input-file stats-dada2-$SOURCE.qza \
  --o-visualization stats-dada2-$SOURCE.qzv


# assign taxonomy 
  qiime feature-classifier classify-sklearn \
    --i-reads rep-seqs-$SOURCE.qza \
    --i-classifier $DB \
    --p-n-jobs 6 \
    --o-classification taxonomy-$SOURCE.qza


# summerize taxonomy info
qiime metadata tabulate \
  --m-input-file taxonomy-$SOURCE.qza \
  --o-visualization taxonomy-$SOURCE.qzv


## remove non-bacterial sequences 

  qiime taxa filter-table \
  --i-table table-$SOURCE.qza \
  --i-taxonomy taxonomy-$SOURCE.qza \
  --p-include Bacteria \
  --p-exclude archaea,eukaryota,mitochondria,chloroplast \
  --o-filtered-table table-bacteria-only-$SOURCE.qza


  qiime taxa filter-seqs \
  --i-sequences rep-seqs-$SOURCE.qza \
  --i-taxonomy taxonomy-$SOURCE.qza \
  --p-include Bacteria \
  --p-exclude archaea,eukaryota,mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-bacteria-only-$SOURCE.qza


## remove singletons <- there should be none as DADA2 removes them!

  qiime feature-table filter-features \
  --i-table table-bacteria-only-$SOURCE.qza \
  --p-min-samples 2 \
  --o-filtered-table table-filtered-$SOURCE.qza


# FeatureTable and FeatureData summarize
qiime feature-table summarize \
  --i-table table-filtered-$SOURCE.qza \
  --o-visualization table-filtered-$SOURCE.qzv \
  --m-sample-metadata-file metadata-$SOURCE.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-bacteria-only-$SOURCE.qza \
  --o-visualization rep-seqs-bacteria-only-$SOURCE.qzv


qiime tools export \
  --input-path table-filtered-$SOURCE.qza \
  --output-path exported-table

biom convert -i exported-table/feature-table.biom -o feature-table-filtered.tsv --to-tsv



### rarefy sequences

# check the number of seqs per sample and set a sampling-depth that the total frequency that each sample should be rarefied to
D=6600

   qiime feature-table rarefy \
    --i-table table-filtered-$SOURCE.qza \
    --p-sampling-depth $D \
    --o-rarefied-table table-rarefied-$SOURCE.qza
    
# export table
qiime tools export \
  --input-path table-rarefied-$SOURCE.qza \
  --output-path ./

biom convert -i feature-table.biom -o feature-table-rarefied-$SOURCE.tsv --to-tsv


### Make a phylogenetic tree
  
 qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 12 \
  --i-sequences rep-seqs-bacteria-only-$SOURCE.qza \
  --o-alignment aligned-rep-seqs-$SOURCE.qza \
  --o-masked-alignment masked-aligned-rep-seqs-$SOURCE.qza \
  --o-tree tree-unrooted-$SOURCE.qza \
  --o-rooted-tree tree-rooted-$SOURCE.qza
  
 qiime tools export \
    --input-path tree-rooted-$SOURCE.qza \
    --output-path ./
    
 mv tree.nwk tree-rooted-$SOURCE.nwk
 

### diversity analysis
# diversity analysis:  # check --p-sampling-depth in the table-filtered.qzv
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny tree-rooted-$SOURCE.qza \
  --i-table table-rarefied-$SOURCE.qza \
  --p-sampling-depth $D \
  --m-metadata-file metadata-$SOURCE.tsv \
  --output-dir core-metrics-results

# exporting matrix
qiime tools export \
  --input-path core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --output-path exported-matrix

mv exported-matrix/distance-matrix.tsv ./unweighted_unifrac_distance_matrix.tsv

qiime tools export \
  --input-path core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --output-path exported-matrix

mv exported-matrix/distance-matrix.tsv ./weighted_unifrac_distance_matrix.tsv

# generate a rarefaction curve
qiime diversity alpha-rarefaction \
  --i-table table-filtered-$SOURCE.qza \
  --i-phylogeny tree-rooted-$SOURCE.qza \
  --p-max-depth 20000 \
  --m-metadata-file metadata-$SOURCE.tsv \
  --o-visualization rarefaction-$SOURCE.qzv \
  --verbose





# well done!
