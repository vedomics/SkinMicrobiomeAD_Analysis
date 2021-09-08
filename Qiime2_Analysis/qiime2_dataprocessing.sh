source activate qiime2-2019.1 # qii version used: 2019-01
cd /path/to/working/directory

## import data
# before running mkdir data_links and add links to each R1 file in casava format. see here: https://docs.qiime2.org/2021.4/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path data_links \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path data.demux_se.qza

# check what we got
qiime demux summarize \
  --i-data data.demux_se.qza \
  --o-visualization data.demux_se.qzv

## remove primer (use only SE data here, thus need to remove only one)
# remove biological primer (27F...leading all seq!)
qiime cutadapt trim-single \
 --p-cores 8 \
 --i-demultiplexed-sequences data.demux_se.qza \
 --p-front AGAGTTTGATCMTGGCTCAG \
 --o-trimmed-sequences data.demux_se.trim.qza

# check what we got
qiime demux summarize \
  --i-data data.demux_se.trim.qza \
  --o-visualization data.demux_se.trim.qzv


##################
###### qc and denoise using deblur
##################

# remove low seq quality reads
qiime quality-filter q-score \
 --i-demux data.demux_se.trim.qza \
 --p-min-quality 4 \
 --o-filtered-sequences data.demux_se.trim.sq4.qza \
 --o-filter-stats filter_stats_sq4.qza

 # check what we got
 qiime demux summarize \
   --i-data data.demux_se.trim.sq4.qza \
   --o-visualization data.demux_se.trim.sq4.qzv

# NOTE 8 cores
# trim 209 for high qual reads only
qiime deblur denoise-16S \
  --i-demultiplexed-seqs data.demux_se.trim.sq4.qza \
  --p-trim-length 209 \
  --p-jobs-to-start 8 \
  --o-representative-sequences data.demux_se.trim.sq4.deblur.qza \
  --o-table deblur-table.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza

# visualize results 
qiime metadata tabulate \
  --m-input-file filter_stats_sq4.qza \
  --o-visualization filter_stats_sq4.qzv
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv

# FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table deblur-table.qza \
  --o-visualization deblur-table.qzv \
  --m-sample-metadata-file metatable_sra_16s_public.tsv 

qiime feature-table tabulate-seqs \
  --i-data data.demux_se.trim.sq4.deblur.qza \
  --o-visualization data.demux_se.trim.sq4.deblur.qzv


########
#### classification using custom-made classifier
########
# classifier build described in paper and available on the lieberman lab website
echo "START classification"
qiime feature-classifier classify-sklearn \
  --i-classifier classifier_27F-534R_100pASV_silva132_min100_max700_trunc209_noMislabelStaph_StaphSpecLvl_rmAmbigous.qza \
  --i-reads data.demux_se.trim.sq4.deblur.qza \
  --o-classification taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qza \
  --o-visualization taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qzv

## barplots
qiime taxa barplot \
  --i-table deblur-table.qza \
  --i-taxonomy taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qza \
  --m-metadata-file metatable_sra_16s_public.tsv \
  --o-visualization taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.bar-plots.qzv

echo "END classification"

# export tax data for analysis
qiime tools export \
  --input-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.bar-plots.qzv \
  --output-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.bar-plots_export

echo "Tax assignment data exported!!!"

##############################################################################
############# Generate a tree for phylogenetic diversity analyses
############# Alpha and beta diversity analysis
##############################################################################

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences data.demux_se.trim.sq4.deblur.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table deblur-table.qza \
  --p-sampling-depth 100 \
  --m-metadata-file metatable_sra_16s_public.tsv \
  --output-dir core-metrics-results_splDepth100

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table deblur-table.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file metatable_sra_16s_public.tsv \
  --output-dir core-metrics-results_splDepth1000


## Explore the microbial composition of the samples in the context of the sample metadata
# 100 reads min
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_splDepth100/faith_pd_vector.qza \
  --m-metadata-file metatable_sra_16s_public.tsv \
  --o-visualization core-metrics-results_splDepth100/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_splDepth100/evenness_vector.qza \
  --m-metadata-file metatable_sra_16s_public.tsv \
  --o-visualization core-metrics-results_splDepth100/evenness-group-significance.qzv

# 1000 reads min
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_splDepth1000/faith_pd_vector.qza \
  --m-metadata-file metatable_sra_16s_public.tsv \
  --o-visualization core-metrics-results_splDepth1000/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_splDepth1000/evenness_vector.qza \
  --m-metadata-file metatable_sra_16s_public.tsv \
  --o-visualization core-metrics-results_splDepth1000/evenness-group-significance.qzv

# export bray curtis dissim matrix for further analysis
qiime tools export \
  --input-path core-metrics-results_splDepth1000/bray_curtis_distance_matrix.qza \
  --output-path core-metrics-results_splDepth1000/bray_curtis_distance_matrix

