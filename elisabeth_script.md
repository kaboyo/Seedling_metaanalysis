# Procedure for processing amplicon sequence data

# Access your hpc 
```
ssh username@iridis5.soton.ac.uk 
```

# Install qiime 2: https://docs.qiime2.org/2024.10/install/index.html
```
# Native installation: 
conda env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
```
# Import the reads (qiime2)

# Prepare the manifest file and transfer files to the cluster

# Download the reads 
https://github.com/rpetit3/fastq-dl

# perform fastp on all the reads (fastp) 
https://www.biorxiv.org/content/10.1101/274100v1.full.pdf)
# Tutorial: 
https://github.com/OpenGene/fastp

# command for fastp
python /scratch/eo2r24/SEEDLING_DATA/data/parallel.py -i /scratch/eo2r24/SEEDLING_DATA/data -o /scratch/eo2r24/SEEDLING_DATA/data/data-output -r /scratch/eo2r24/SEEDLING_DATA/data/reports -a '-f 3 -t 2'

# Transfer the file from the cluster to the local computer
scp eo2r24@iridis5.soton.ac.uk:/scratch/eo2r24/SEEDLING_DATA/data/current-files.txt /Users/expeditoolimi/Documents/Southampton/Meta-analysis

# Transfer the file from the local to cluster
scp import_manifest_cluster eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/import_manifest_cluster

# Download the metadata
fastq-dl --accession  PRJNA667691 --only-download-metadata 

2. Importing data in qiime2 (using the QIIME2 commands directly on the terminal)
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path import_manifest_cluster --input-format PairedEndFastqManifestPhred33 --output-path /home/eo2r24/QIIME2_TomaSeed_Olimi/TomatoSeed1000_demux.qza
#######################
  OR
####For loop command to import the files
for file in *.fq.gz; do
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path import_manifest_cluster \
      --output-path ~/QIIME2_TomaSeed_Olimi/TomatoSeed1000_demux.qza \
      --input-format PairedEndFastqManifestPhred33
done

command: 
for file in *.fq.gz; do     qiime tools import       --type 'SampleData[PairedEndSequencesWithQuality]' --input-path import_manifest_cluster --output-path ~/QIIME2_TomaSeed_Olimi/TomatoSeed1000_demux.qza --input-format PairedEndFastqManifestPhred33; done

###############################################
#### Quality Check
#On the local computer
qiime demux summarize --i-data /Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/TomatoSeed1000_demux.qza --o-visualization /Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/TomatoSeed1000_demux.qzv


####On the cluster
qiime demux summarize --i-data TomatoSeed1000_demux.qza --o-visualization /home/eo2r24/QIIME2_TomaSeed_Olimi/TomatoSeed1000_demux.qzv

qiime demux summarize --i-data TomatoSeed1000_demux.qza --o-visualization TomatoSeed1000_demux.qzv 


#Visualize the output (# check for Qiime2 view for quality and verify length.)
qiime tools view /Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/TomatoSeed1000_demux.qzv


################################################################################################################################################################################
Note that first dada had no timing or truncation

# 3. Run DADA2 dor demoising the reads
###Option 1 V4 Primers Earth Microbiome project

####Processed on R version 4.2.2 (2022-10-31) 
####DADA2: 1.26.0 / Rcpp: 1.0.11 / RcppParallel: 5.1.6 
#####package ‘optparse’ was built under R version 4.2.3 
#####Loading required package: Rcpp

##no trim pr truncation of the reads
Command on the (computer)#1
qiime dada2 denoise-paired --i-demultiplexed-seqs /Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/TomatoSeed1000_demux.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads 2 --o-representative-sequences /Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/DADA2/TomatoSeed_repSeq_1.qza --o-table /Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/DADA2/TomatoSeed_table_1.qza --o-denoising-stats /Users/expedxitoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/DADA2/denoiseStats_1.qza --verbose


##--p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 180 --p-trunc-len-r 180

qiime dada2 denoise-paired --i-demultiplexed-seqs TomatoSeed1000_demux.qza --p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 180 --p-trunc-len-r 180 --p-n-threads 3 --o-representative-sequences DADA2/TomatoSeed_repSeq180.qza --o-table DADA2/TomatoSeed_table180.qza --o-denoising-stats DADA2/denoiseStats180.qza --verbose

Export to remote
scp /Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Olimi_analysis_Global_seed_microbiome/QIIME/DADA2 eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2


scp denoiseStats_1.qzv eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/Notrunctrim
scp denoiseStats_1.qza eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/Notrunctrim
scp TomatoSeed_table_1.qza eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/Notrunctrim
scp TomatoSeed_repSeq_1.qza eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/Notrunctrim

scp Veren.sh eo2r24@iridis5.soton.ac.uk:~/scratch/eo2r24/Verena-data
#######################################################################################
#On the cluster

Warning message:
package 'optparse' was built under R version 4.2.3 
R version 4.2.2 (2022-10-31) 
Loading required package: Rcpp
DADA2: 1.26.0 / Rcpp: 1.0.11 / RcppParallel: 5.1.6 


Command  (on cluster:sinteractive --partition=highmem --cpus-per-task=30 --mem=600000 --nodes=5 --time=01-02:03:04)#2

##--p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads
                 
qiime dada2 denoise-paired --i-demultiplexed-seqs TomatoSeed1000_demux.qza --p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads 40 --o-representative-sequences DADA2/TomatoSeed_repSeq.qza --o-table DADA2/TomatoSeed_table.qza --o-denoising-stats DADA2/denoiseStats.qza --verbose

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
https://github.com/benjjneb/dada2/issues/213:cutadapt
#enter primer sequences
fwdp <-"GTGYCAGCMGCCGCGGTAA" #forward
revp <-"GGACTACNVGGGTWTCTAAT" #reverse

cat("#!/bin/bash \n",file="cutadapt.sh")
cat("\nexport _JAVA_OPTIONS=\"-Xmx16G\"\n",file="cutadapt.sh",append=TRUE)
cat("\nfor FWD in $(ls *R1*); do",file="cutadapt.sh",append=TRUE)
cat("\n\tBASE=$(echo $FWD | cut -d \"_\" -f1);",file="cutadapt.sh",append=TRUE)
cat("\n\tREV=$(echo $FWD | sed 's/R1/R2/g');",file="cutadapt.sh",append=TRUE)
cat("\n\tcutadapt -g ", fwdp, " -G ",revp, " -o $BASE\"_P1.fastq.gz\" -p $BASE\"_P2.fastq.gz\" $FWD $REV;",file="cutadapt.sh",append=TRUE)
cat("\ndone",file="cutadapt.sh",append=TRUE)
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################

##command to visualise the statistics:  
qiime metadata tabulate --m-input-file denoiseStats.qza --o-visualization denoiseStats.qzv
qiime tools view denoiseStats.qzv

###############################################################################################################################################################################
Performed on the cluster with parameters: sinteractive --partition=highmem --cpus-per-task=30 --mem=600000 --nodes=5 --time=01-02:03:04

qiime feature-classifier classify-consensus-vsearch --i-query TomatoSeed_repSeq180.qza --i-reference-reads silva-138-99-seqs.qza --i-reference-taxonomy silva-138-99-tax.qza --p-threads 40 --o-classification TomatoSeed_taxonomy.qza --o-search-results TomatoSeed_tophits.qza --output-dir Taxonomy/ --verbose

#view the feature table and taxonomony
###Feature table
generate: qiime feature-table summarize --i-table TomatoSeed_table180.qza --o-visualization TomatoSeed_table180.qzv 
Import from the remote cluster:#scp eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/TomatoSeed_table180.qzv ~/quality
visualize: qiime tools view TomatoSeed_table180.qzv

###Taxonomy
generate a visualization of taxonomy: qiime metadata tabulate --m-input-file TomatoSeed_taxonomy.qza --o-visualization TomatoSeed_taxonomy.qzv
Import from the remote cluster:#scp eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/TomatoSeed_taxonomy.qzv ~/quality
visualize: qiime tools view TomatoSeed_taxonomy.qzv


##filter chloro and mitochondria and unassigned phlum

qiime taxa filter-table --i-table TomatoSeed_table180.qza --i-taxonomy TomatoSeed_taxonomy.qza --p-include p__  --p-exclude chloroplast,mitochondria --o-filtered-table TomatoSeed_table_noMito_noChloro.qza --verbose

Remove archaea

qiime taxa filter-table --i-table TomatoSeed_table_noMito_noChloro.qza --i-taxonomy TomatoSeed_taxonomy.qza --p-exclude “d__Archaea” --o-filtered-table TomatoSeed_table_noMito_noChloro_archaea.qza --verbose


#Visualize the filtered tables
generate:

qiime feature-table summarize --i-table TomatoSeed_table_noMito_noChloro.qza --o-visualization TomatoSeed_table_noMito_noChloro.qza.qzv


qiime feature-table summarize --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --o-visualization TomatoSeed_table_noMito_noChloro_archaea.qzv


Import from remote: 
scp eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/TomatoSeed_table_noMito_noChloro.qzv ~/quality/DADA2

scp eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2 TomatoSeed_table_noMito_noChloro_archaea.qzv ~/quality/DADA2
scp eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2 TomatoSeed_table_noMito_noChloro_archaea.qza ~/quality/DADA2
View: 
qiime tools view TomatoSeed_table_noMito_noChloro.qzv

qiime tools view TomatoSeed_table_noMito_noChloro_archaea.qzv

##################################################################################################################################
#Obtain the different rarefactions the different maximum depths
scp Metadata_rarefaction.tsv eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2

#1500: 
qiime diversity alpha-rarefaction --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-metadata-file Metadata_rarefaction.tsv  --p-max-depth 1500  --output-dir  rarefy_1500/  --verbose
qiime tools view visualization.qzv 

#2000: 
qiime diversity alpha-rarefaction --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-metadata-file Metadata_rarefaction.tsv  --p-max-depth 2000  --output-dir  rarefy_2000/  --verbose
qiime tools view visualization.qzv 

#22000: 
qiime diversity alpha-rarefaction --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-metadata-file Metadata_rarefaction.tsv  --p-max-depth 22000  --output-dir  rarefy_22000/  --verbose
qiime tools view visualization.qzv 
#220000: 
qiime diversity alpha-rarefaction --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-metadata-file Metadata_rarefaction.tsv  --p-max-depth 220000  --output-dir  rarefy_220000/  --verbose
qiime tools view visualization.qzv 
#8000: 
qiime diversity alpha-rarefaction --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-metadata-file Metadata_rarefaction.tsv  --p-max-depth 8000  --output-dir  rarefy_8000/  --verbose
qiime tools view visualization.qzv 

Intuition: 
Generate interactive alpha rarefaction curves by computing rarefactions between `min_depth` and `max_depth`. The number of intermediate depths to compute is controlled by the `steps` parameter, with n `iterations` being computed at each rarefaction depth. If sample metadata is provided, samples may be grouped based on distinct values within a metadata column.
##################################################################################################################################

##core500 (no phylogenetic)
qiime diversity core-metrics --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --p-sampling-depth 500 --m-metadata-file Metadata.tsv --output-dir core/

qiime feature-table summarize --i-table rarefied_table.qza --o-visualization rarefied_table.qzv 
qiime tools view rarefied_table.qzv

##70 samples verloren-..
qiime taxa barplot --i-table rarefied_table.qza --i-taxonomy /Users/expeditoolimi/quality/DADA2/TomatoSeed_taxonomy.qza  --m-metadata-file Metadata.tsv --o-visualization Barplot_rarefied_table.qzv

qiime tools view Barplot_rarefied_table.qzv


###Filter the representative sequences

qiime feature-table filter-seqs --i-data TomatoSeed_repSeq180.qza --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --o-filtered-data TomatoSeed_filtered_repSeq180.qza --verbose 


#https://zenodo.org/records/10531254 #coverM


##tree without archaea

Export the TomatoSeed_repSeq-noMito_noChloro.qza into remote disk
scp TomatoSeed_table_noMito_noChloro_archaea.qza eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2
scp TomatoSeed_filtered_repSeq180.qza eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2


qiime phylogeny align-to-tree-mafft-fasttree --i-sequences TomatoSeed_filtered_repSeq180.qza --output-dir tree/ --p-n-threads 50 --verbose

#Tree files
Saved FeatureData[AlignedSequence] to: tree/alignment.qza
Saved FeatureData[AlignedSequence] to: tree/masked_alignment.qza
Saved Phylogeny[Unrooted] to: tree/tree.qza
Saved Phylogeny[Rooted] to: tree/rooted_tree.qza

Import the tree files
scp -r eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/tree ~/quality/DADA2

###########################################################################################################################################
https://holtlab.net/2015/10/03/plotting-trees-data/: the tree plotting software
###########################################################################################################################################
########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
##################################################################################################
Heatmap
Weighted (—p-method)
qiime feature-table heatmap --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-sample-metadata-file Metadata.tsv --m-sample-metadata-column Cultivar_summary --p-title Tree-of-cultivars --p-metric braycurtis --p-method weighted --p-color-scheme coolwarm --output-dir heatmap_weighted/ --verbose

Weighted (—p-method)
qiime feature-table heatmap --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-sample-metadata-file Metadata.tsv --m-sample-metadata-column Cultivar_summary --p-title Tree-of-cultivars --p-metric braycurtis --p-method average --p-color-scheme coolwarm --output-dir heatmap_average/ --verbose

Weighted (—p-method)
qiime feature-table heatmap --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-sample-metadata-file Metadata.tsv --m-sample-metadata-column Cultivar_summary --p-title Tree-of-cultivars --p-metric braycurtis --p-method centroid --p-color-scheme coolwarm --output-dir heatmap_centroid/ --verbose

Weighted (—p-method)
qiime feature-table heatmap --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-sample-metadata-file Metadata.tsv --m-sample-metadata-column Cultivar_summary --p-title Tree-of-cultivars --p-metric braycurtis --p-method complete --p-color-scheme coolwarm --output-dir heatmap_complete/ --verbose

Weighted (—p-method)
qiime feature-table heatmap --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --m-sample-metadata-file Metadata.tsv --m-sample-metadata-column Cultivar_summary --p-title Tree-of-cultivars --p-metric braycurtis --p-method ward --p-color-scheme coolwarm --output-dir heatmap_ward/ --verbose

#########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#################################################################################################

Transfer whole tree folder on to the computer
Command: scp -r eo2r24@iridis5.soton.ac.uk:~/QIIME2_TomaSeed_Olimi/DADA2/tree ~/quality/DADA2

##core phylogeny

qiime diversity core-metrics-phylogenetic --i-table TomatoSeed_table_noMito_noChloro_archaea.qza --i-phylogeny tree/rooted_tree.qza --p-sampling-depth 500 --m-metadata-file Metadata.tsv --output-dir core500_phylogeny/

-qiime tools view bray_curtis_emperor.qzv
-qiime tools view jaccard_emperor.qzv
-qiime tools view unweighted_unifrac_emperor.qzv
-qiime tools view weighted_unifrac_emperor.qzv
##################################################################################################
##################################################################################################
##################################################################################################
qiime taxa barplot --i-table rarefied_table.qza --i-taxonomy /Users/expeditoolimi/quality/DADA2/TomatoSeed_taxonomy.qza --m-metadata-file /Users/expeditoolimi/quality/DADA2/Metadata.tsv --o-visualization Barplot_rarefied_table.qzv

qiime tools view Barplot_rarefied_table.qzv
#######################
##grouping
qiime feature-table group --i-table rarefied_table.qza --p-axis sample --m-metadata-column English_Name --p-mode mean-ceiling --o-grouped-table rarefied_table_groupedCultivar.qza --m-metadata-file /Users/expeditoolimi/quality/DADA2/Metadata.tsv


qiime feature-table summarize --i-table rarefied_table_groupedCultivar.qza --o-visualization rarefied_table_groupedCultivar.qzv 

qiime tools view rarefied_table_groupedCultivar.qzv 

#######################

qiime taxa barplot --i-table rarefied_table_groupedCultivar.qza --i-taxonomy /Users/expeditoolimi/quality/DADA2/TomatoSeed_taxonomy.qza --o-visualization Barplot_grouped_rarefied_table.qzv --m-metadata-file Metadata.tsv
ERROR: failed..
##################################################################################################

##stats
#shannon
qiime diversity alpha-group-significance --i-alpha-diversity shannon_vector.qza  --m-metadata-file Metadata.tsv --o-visualization alpha-group-sig_shannon.qzv

qiime tools view  alpha-group-sig_shannon.qzv 


#pd 

qiime diversity alpha-group-significance --i-alpha-diversity faith_pd_vector.qza  --m-metadata-file Metadata.tsv --o-visualization alpha-group-sig_faithpd.qzv

qiime tools view  alpha-group-sig_faithpd.qzv 


#Observed ASVs
qiime diversity alpha-group-significance --i-alpha-diversity  observed_features_vector.qza  --m-metadata-file Metadata.tsv --o-visualization alpha-group-sig_observed_features.qzv

qiime tools view  alpha-group-sig_observed_features.qzv 

evenness
qiime diversity alpha-group-significance --i-alpha-diversity  evenness_vector.qza  --m-metadata-file Metadata.tsv --o-visualization alpha-group-sig_evenness_features.qzv

qiime tools view  alpha-group-sig_evenness_features.qzv 

##################################################################################################
##################################################################################################

Adonis 

#Bray-curtis
qiime diversity adonis --i-distance-matrix  bray_curtis_distance_matrix.qza --m-metadata-file Metadata.tsv --p-formula Cultivar_summary*Country_of_Origin*Production_site*Fruit_shape*Fruit_color*Fruit_taste --o-visualization adonis_BrayCurtis.qzv

qiime tools view  adonis_BrayCurtis.qzv 

######################
#Jaccard
qiime diversity adonis --i-distance-matrix  jaccard_distance_matrix.qza --m-metadata-file Metadata.tsv --p-formula Cultivar_summary*Country_of_Origin*Production_site*Fruit_shape*Fruit_color*Fruit_taste --o-visualization adonis_jaccard.qzv

qiime tools view  adonis_jaccard.qzv

######################
#unweighted_unifrac_distance_matrix.qza

qiime diversity adonis --i-distance-matrix  unweighted_unifrac_distance_matrix.qza --m-metadata-file Metadata.tsv --p-formula Cultivar_summary*Country_of_Origin*Production_site*Fruit_shape*Fruit_color*Fruit_taste --o-visualization adonis_unweighted_unifrac_distance.qzv

qiime tools view  adonis_adonis_unweighted_unifrac_distance.qzv

######################
#weighted_unifrac_distance_matrix.qza
qiime diversity adonis --i-distance-matrix  weighted_unifrac_distance_matrix.qza --m-metadata-file Metadata.tsv --p-formula Cultivar_summary*Country_of_Origin*Production_site*Fruit_shape*Fruit_color*Fruit_taste --o-visualization adonis_weighted_unifrac_distance.qzv

qiime tools view  adonis_adonis_weighted_unifrac_distance.qzv

######################
#weighted_unifrac_distance_matrix.qza
qiime diversity adonis --i-distance-matrix  weighted_unifrac_distance_matrix.qza --m-metadata-file Metadata.tsv --p-formula Cultivar_summary*Country_of_Origin*Production_site*Fruit_shape*Fruit_color*Fruit_taste --o-visualization adonis_weighted_unifrac_distance.qzv
qiime tools view  adonis_adonis_weighted_unifrac_distance.qzv

##################################################################################################
#Plotting geocoordinates colored by alpha diversity values
#To be done in R
##################################################################################################
##core
#Default (0.5)

qiime feature-table core-features --i-table  rarefied_table.qza --p-min-fraction 0.5 --o-visualization Core_features_rarefied_table05.qzv

qiime tools view Core_features_rarefied_table05.qzv 

#0.1
qiime feature-table core-features --i-table  rarefied_table.qza --p-min-fraction 0.1 --o-visualization Core_features_rarefied_table01.qzv

qiime tools view Core_features_rarefied_table01.qzv


##biplot

qiime feature-table relative-frequency --i-table rarefied_table.qza --o-relative-frequency-table table_relative_frequencies.qza

qiime feature-table summarize --i-table table_relative_frequencies.qza --o-visualization table_relative_frequencies.qzv


qiime diversity pcoa-biplot --i-pcoa weighted_unifrac_pcoa_results.qza --i-features table_relative_frequencies.qza --o-biplot biplot_WUnif.qza

qiime diversity pcoa-biplot --i-pcoa bray_curtis_pcoa_results.qza --i-features table_relative_frequencies.qza --o-biplot biplot_bray.qza


qiime emperor biplot --i-biplot biplot_WUnif.qza --m-sample-metadata-file Metadata.tsv --o-visualization biplot_WUnif.qzv
qiime tools view biplot_WUnif.qzv

qiime emperor biplot --i-biplot biplot_bray.qza --m-sample-metadata-file Metadata.tsv --o-visualization biplot_bray.qzv
qiime tools view biplot_bray.qzv


##export

qiime tools export --input-path rarefied_table.qza --output-path exported_table_rarefied_500depth

biom convert -i feature-table.biom --to-tsv -o feature-table_filtered_rarefied500depth.txt

qiime tools export --input-path TomatoSeed_taxonomy.qza --output-path Taxonomy_exported

qiime tools export --input-path tree/tree.qza --output-path Taxonomy_exported/

Export the filtered table

qiime tools export --input-path TomatoSeed_table_noMito_noChloro_archaea.qza --output-path Taxonomy_exported/

biom convert -i feature-table.biom --to-tsv -o feature-table_filtered_rarefied500depth.txt
biom convert -i feature-table.biom --to-tsv -o feature-table_filtered.txt

#########################################################################################
#########################################################################################
#########################################################################################
###table group -- rel frequency

qiime feature-table group --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2/filteredUnassigned/core500_phylogeny/rarefied_table.qza' --p-axis sample --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/mappingFile_TomatoSeed100' --m-metadata-column ProductionSite --p-mode mean-ceiling --o-grouped-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2/filteredUnassigned/core500_phylogeny/table_group/Table_grouped_ProductionSite.qza'

qiime feature-table relative-frequency --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2/filteredUnassigned/core500_phylogeny/table_group/Table_grouped_ProductionSite.qza' --o-relative-frequency-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2/filteredUnassigned/core500_phylogeny/table_group/Table_grouped_ProductionSite_relFreq.qza'
#########################################################################################
#########################################################################################
#########################################################################################

###import Simonin et al., 2021 data

biom convert -i '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S_table-FINAL-V4-MiSeq-filtered.txt' -o '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S_table-FINAL-V4-MiSeq-filtered.biom' --table-type="OTU table" --to-json


qiime tools import --input-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S_table-FINAL-V4-MiSeq-filtered.biom' --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S_table-FINAL-V4-MiSeq-filtered.qza'

qiime tools import --input-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-filtered.fasta' --output-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-filtered.qza' --type 'FeatureData[Sequence]'

##filter tomato studies
qiime feature-table filter-samples --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S_table-FINAL-V4-MiSeq-filtered.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/MappingFile_all_modified' --p-where "Plant='Tomato'" --o-filtered-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/50Seed_Tomato_16S.qza' 

core_tomato part barret
qiime diversity core-metrics --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/50Seed_Tomato_16S.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/mappingfile_tomato_part' --p-sampling-depth 500 --output-dir '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/core500_tomato'


##merge with tomato100
qiime feature-table merge --i-tables '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/50Seed_Tomato_16S.qza' --i-tables '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2/TomatoSeed_table_noMito_noChloro_noUnassigned.qza' --o-merged-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/Table_Tomato_all.qza'

qiime feature-table merge-seqs --i-data '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2/TomatoSeed_repSeq-noMito_noChloro.qza' --i-data '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-filtered.qza' --o-merged-data '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/RepSeqs_Tomato_all.qza'

##core joined with tomatoseed100
qiime diversity core-metrics --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/Table_Tomato_all.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/mappingFile_joined' --p-sampling-depth 500 --output-dir '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined'

##taxonomy together

qiime feature-classifier classify-consensus-vsearch --i-query '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/RepSeqs_Tomato_all.qza' --i-reference-reads '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Databases/vsearch_16S_DB_Silva/rep-set_SILVA132_consensus_99.qza' --i-reference-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Databases/vsearch_16S_DB_Silva/ref-taxonomy_SILVA132_consensus_99.qza' --p-threads 30 --o-classification '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/Taxonomy_silva132_Tomato_all.qza'


##barplot
qiime taxa barplot --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/Table_Tomato_all.qza' --i-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/Taxonomy_silva132_Tomato_all.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/mappingFile_joined' --o-visualization '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/barplot_Tomato_all.qzv' 


qiime taxa barplot --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/rarefied_table.qza' --i-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/Taxonomy_silva132_Tomato_all.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/mappingFile_joined' --o-visualization '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/Barplot_rarefied_table.qzv' 



##group and relative abundance

qiime feature-table group --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/rarefied_table.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/mappingFile_joined' --m-metadata-column Genotype --p-axis sample --p-mode mean-ceiling --o-grouped-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/Table_Genotype_Tomato_all.qza' 


qiime feature-table relative-frequency --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/Table_Genotype_Tomato_all.qza' --o-relative-frequency-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/Table_relFreq_Genotype_Tomato_all.qza'

##export
qiime tools export --input-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/Table_relFreq_Genotype_Tomato_all.qza' --output-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/exported'

biome convert -i '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/exported/feature-table_relFreq_StudyID_Tomato_all.biom' -o '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/exported/feature-table_relFreq_StudyID_Tomato_all.txt' --to-tsv 



##core features
qiime feature-table core-features --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/rarefied_table.qza' --o-visualization '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/Core_features_rarefied_table.qzv' --p-min-fraction 0.1

##collapse lvl6
qiime taxa collapse --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/Table_Genotype_Tomato_all.qza' --i-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/Taxonomy_silva132_Tomato_all.qza' --p-level 6 --o-collapsed-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/core500_joined/grouped_table/Table_relFreq_Genotype_Tomato_all_lvl6.qza'



##because grouping in pcoa is weird (e.g. together with simonin samples), try to cut different in dada2 (according to simonin documentation)
###dada2 with cut 180 
qiime dada2 denoise-paired --i-demultiplexed-seqs '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/TomatoSeed1000_demux.qza' --p-trunc-len-f 180 --p-trunc-len-r 180 --p-trim-left-f 20 --p-trim-left-r 20 --p-n-threads 30 --o-representative-sequences '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_repSeq.qza' --o-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_table.qza' --o-denoising-stats '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/denoiseStats.qza' --verbose

qiime metadata tabulate --m-input-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/denoiseStats.qza' --o-visualization '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/denoiseStats.qzv'


qiime feature-classifier classify-consensus-vsearch --i-query '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_repSeq.qza' --i-reference-reads '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Databases/vsearch_16S_DB_Silva/rep-set_SILVA132_consensus_99.qza' --i-reference-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Databases/vsearch_16S_DB_Silva/ref-taxonomy_SILVA132_consensus_99.qza' --p-threads 30 --o-classification '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_taxonomy_silva132.qza'


qiime taxa filter-table --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_table.qza' --i-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_taxonomy_silva132.qza' --p-exclude chloroplast,mitochondria,unassigned --o-filtered-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_table_noMito_noChloro_noUnass.qza'

qiime feature-table merge --i-tables '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_table_noMito_noChloro_noUnass.qza' --i-tables '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/50Seed_Tomato_16S.qza' --o-merged-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Table_merged_all_tomato.qza' 

qiime diversity core-metrics --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Table_merged_all_tomato.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/merged_with_tomatoseed100/mappingFile_joined' --p-sampling-depth 500 --output-dir '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500' 

##filtering is still the same, but table of new one is better then old...


##filter RepSeqs
qiime taxa filter-seqs --i-sequences '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_repSeq.qza' --i-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_taxonomy_silva132.qza' --p-exclude chloroplast,mitochondria,unassigned --o-filtered-sequences '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_repSeq_noMito_noChloro_noUnass.qza'

##join with simonin
qiime feature-table merge-seqs --i-data '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/TomatoSeed_repSeq_noMito_noChloro_noUnass.qza' --i-data '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/50Seed_data/Subset1-All_studies_merged_16S-rep-seqs-FINAL-V4-MiSeq-filtered.qza' --o-merged-data '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/RepSeqs_filtered_merged_all_tomato.qza'

##tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/RepSeqs_filtered_merged_all_tomato.qza' --output-dir '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/tree_all' --p-n-threads 30 --verbose



##taxonomy
qiime feature-classifier classify-consensus-vsearch --i-query '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/RepSeqs_filtered_merged_all_tomato.qza' --i-reference-reads '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Databases/vsearch_16S_DB_Silva/rep-set_SILVA132_consensus_99.qza' --i-reference-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Databases/vsearch_16S_DB_Silva/ref-taxonomy_SILVA132_consensus_99.qza' --p-threads 30 --o-classification '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Taxonomy_Silva132_merged_all_tomato.qza'


##core-phylogeny
qiime diversity core-metrics-phylogenetic --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Table_merged_all_tomato.qza' --i-phylogeny '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/tree_all/rooted_tree.qza' --p-sampling-depth 500 --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/mappingFile_joined' --output-dir '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic' 



#iTOL

qiime feature-table group --p-axis sample --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/mappingFile_joined' --m-metadata-column Study_ID --p-mode mean-ceiling --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Table_merged_all_tomato.qza' --o-grouped-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Table_grouped_Study_merged_all_tomato.qza' 

qiime tools export --input-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Table_grouped_Study_merged_all_tomato.qza' --output-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/exported'

biom convert -i '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/exported/feature-table_grouped_StudyID_merged_tomato.biom' -o '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/exported/feature-table_grouped_StudyID_merged_tomato.txt' --to-tsv 

qiime tools export --input-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Taxonomy_Silva132_merged_all_tomato.qza' --output-path '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/exported'



##
qiime taxa collapse --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table.qza' --i-taxonomy '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/Taxonomy_Silva132_merged_all_tomato.qza' --p-level 6 --o-collapsed-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table_lvl6.qza'

qiime feature-table relative-frequency --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table_lvl6.qza' --o-relative-frequency-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table_lvl6_relFreq.qza' 

qiime feature-table group --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table_lvl6.qza' --m-metadata-file '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/mappingFile_joined' --m-metadata-column Study_ID --p-axis sample --p-mode mean-ceiling --o-grouped-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table_lvl6_grouped_Study.qza' 

qiime feature-table relative-frequency --i-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table_lvl6_grouped_Study.qza' --o-relative-frequency-table '/mnt/9cd04016-2a6f-4323-a7f8-e462540180c7/Win10_QIIME_Data_Folder/Peter/Tomi_TomatoSeed1000/qiime2/dada2_cut180/merged_with_Simonin/core500_phylogenetic/rarefied_table_lvl6_grouped_Study_relFreq.qza' 

##export for comparison on core features
