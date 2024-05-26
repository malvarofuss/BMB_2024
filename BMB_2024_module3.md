---
layout: aws_tutorial_page
permalink: /BMB_2024_module3
title: AWS 2024 - Module 3
header1: Workshop Pages for Students
header2: Module 3 - Microbiome statistics and visualizations
image: /site_images/BMB_2024_v1.png
home: https://bioinformaticsdotca.github.io/BMB_2024
---

This tutorial is part of the 2024 Canadian Bioinformatics Workshops Beginner Microbiome Analysis (St John's, NL, May 27-28). It is based on the Amplicon SOP v2 available on the [Microbiome Helper repository](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)) and previous workshops designed by Diana Haider and Robert Beiko.

**Author**: Monica Alvaro Fuss

## Table of Contents

1.  [Build tree](#1-build-tree-withsepp-qiime-2-plugin)
2.  [Generate rarefaction curves](#2-generate-rarefaction-curves)
3.  [Calculating diversity metrics and generating ordination plots](#3-calculating-diversity-metrics-and-generating-ordination-plots)
4.  [Generate stacked bar chart of taxa relative abundances](#4-generate-stacked-bar-chart-of-taxa-relative-abundances)
5.  [Identifying differentially abundant features with ANCOM](#5-identifying-differentially-abundant-features-with-ancom)
6.  [Exporting the final abundance, profile and sequence files](#6-exporting-the-final-abundance-profile-and-sequence-files)

## Introduction

Modules 2 and 3 provide a walkthrough of an end-to-end pipeline using the command line interface for the analysis of high-throughput marker gene data. Commonly used marker genes for microbiome analysis include the 16S ribosomal RNA (rRNA) for prokaryotes, 18S rRNA for eukaryotes, and the internal transcribed spacer (ITS) for fungi. In this tutorial, we will explore a 16S rRNA dataset from wild blueberry (*Vaccinium angustifolium*) soil communities from both bulk and rhizosphere soils associated with either natural or managed habitats. You can read more about the study in these papers:

-   [Variation in Bacterial and Eukaryotic Communities Associated with Natural and Managed Wild Blueberry Habitats](https://apsjournals.apsnet.org/doi/10.1094/PBIOMES-03-17-0012-R)
-   [Metagenomic Functional Shifts to Plant Induced Environmental Changes](https://www.frontiersin.org/articles/10.3389/fmicb.2019.01682/full#B50)

In Module 2 we covered the basics of marker gene analysis from raw reads to filtered feature table, and in Module 3 we will explore examples of downstream analyses that can be used to draw biological conclusions from these data. The pipeline described is embedded in the latest version of QIIME2 (Quantitative Insights into Microbial Ecology version 2023.2), which is a popular microbiome bioinformatics platform for microbial ecology built on user-made software packages called plugins that work on QIIME2 artifact or QZA files. Documentation for these plugins can be found in the [QIIME 2 user documentation](https://docs.qiime2.org/2023.2/), along with tutorials and other useful information. QIIME2 also provides interpretable visualizations that can be accessed by opening any generated QZV files within [QIIME2 View](https://view.qiime2.org/).

Let's set up our Module 3 directory inside `workspace` and link to the data we will be using.

```         
cd ~/workspace
mkdir bmb_module3
cd bmb_module3
ln -s ~/CourseData/MIC_data/BMB_data/deblur_output/ .
ln -s ~/CourseData/MIC_data/BMB_data/taxa .
ln -s ~/CourseData/MIC_data/BMB_data/Blueberry_metadata_reduced.tsv .
```

If you deactivated your QIIME2 environment, reactivate it with either of the commands below.

```         
conda activate qiime2-amplicon-2024.2
```

```         
conda activate qiime2-amplicon-backup
```

When you are finished this tutorial you can deactivate the conda environment using:

```         
conda deactivate
```

Throughout this module, there are some questions aimed to help your understanding of some of the key concepts. You'll find the answers at the bottom of this page, but no one will be marking them.

## 1. Build tree with [SEPP QIIME 2 plugin](https://github.com/qiime2/q2-fragment-insertion)

[SEPP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5904434/) is one method for placing short sequences into a reference phylogenetic tree. This is a useful way of determining a phylogenetic tree for your ASVs. **For 16S data** you can do this with q2-fragment-insertion using the below command.

**Again, due to memory constraints, you can copy the output into your folder with the following command:**

```         
ln -s ~/CourseData/MIC_data/BMB_data/asvs-tree.qza .
ln -s ~/CourseData/MIC_data/BMB_data/insertion-placements.qza .
```

```         
qiime fragment-insertion sepp \
  --i-representative-sequences deblur_output/rep_seqs_final.qza \
  --i-reference-database ~/CourseData/MIC_data/BMB_data/sepp-refs-gg-13-8.qza \
  --o-tree asvs-tree.qza \
  --o-placements insertion-placements.qza \
  --p-threads 4
```

`sepp-refs-gg-13-8.qza` can be downloaded as specified in the [fragment-insertion instructions](https://library.qiime2.org/plugins/q2-fragment-insertion/16/). You can specify custom reference files to place other amplicons, but the easiest approach for **18S and ITS data** is to instead create a *de novo* tree as outlined in the [Microbiome Helper repository](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)).

## 2. Generate rarefaction curves

A key quality control step is to plot rarefaction curves for all of your samples to determine if you performed sufficient sequencing. The below command will generate these plots (make sure you have the correct maximum sequencing depth as per your filtered feature table).

```         
qiime diversity alpha-rarefaction \
  --i-table deblur_output/deblur_table_final.qza \
  --p-max-depth X \
  --p-steps 20 \
  --i-phylogeny asvs-tree.qza \
  --o-visualization rarefaction_curves.qzv
```

**Question 1: What is a good rarefaction depth for diversity analysis?**

## 3. Calculating diversity metrics and generating ordination plots

Common alpha and beta-diversity metrics can be calculated with a single command in QIIME 2. Ordination plots (such as PCoA plots for weighted UniFrac distances) will be generated automatically as well. This command will also rarefy all samples to the sample sequencing depth before calculating these metrics (X is a placeholder for the lowest *reasonable* sample depth; samples with depth below this cut-off will be excluded).

```         
qiime diversity core-metrics-phylogenetic \
  --i-table deblur_output/deblur_table_final.qza \
  --i-phylogeny asvs-tree.qza \
  --p-sampling-depth X  \
  --m-metadata-file Blueberry_metadata_reduced.tsv \
  --p-n-jobs-or-threads 4 \
  --output-dir diversity
```

For alpha diversity visualizations, you will need to produce boxplots comparing the different categories in your metadata file. For example, to create boxplots comparing the Shannon alpha-diversity metric you can use this command:

```         
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity/shannon_vector.qza \
  --m-metadata-file Blueberry_metadata_reduced.tsv \
  --o-visualization diversity/shannon_compare_groups.qzv
```

*Hint: you will need to change this command for the other alpha diversity metrics. You can see the other metrics available by running `ls diversity/*_vector.qza`*.

Note that you can also export (see below) this or any other diversity metric file (ending in *.qza*) and analyze them with a different program.

**Question 2: are there any significant differences in alpha diversity between any of our metadata categories?**

**Question 3: which metadata category appears to provide more separation in the beta diversity PCoA plots?**

If you want to run a PERMANOVA test to calculate differences in beta diversity between groups of interest, you can run the following command. Remember you will need to change it according to which beta diversity metric you are using, and you will also have to specify your metadata category of interest.

```         
qiime diversity beta-group-significance \
  --i-distance-matrix diversity/bray_curtis_distance_matrix.qza \
  --m-metadata-file Blueberry_metadata_reduced.tsv \
  --m-metadata-column Desription_1 \
  --o-visualization diversity/bray_curtis_compare_groups_Description_1.qzv
```

**Question 4: what do you mean I'm getting an error message? It should be working fine, it worked last time! I can't see what I'm doing wro... oooooh there we go, I see now :)**

## 4. Generate stacked bar chart of taxa relative abundances

Another useful output is the interactive stacked bar-charts of the taxonomic abundances across samples, which can be output with this command:

```         
qiime taxa barplot \
  --i-table deblur_output/deblur_table_final.qza \
  --i-taxonomy taxa/classification.qza \
  --m-metadata-file Blueberry_metadata_reduced.tsv \
  --o-visualization taxa/taxa_barplot.qzv
```

**Question 5: can you identify any patterns between the metadata groups?**

## 5. Identifying differentially abundant features with ANCOM

[ANCOM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/) is one method to test for differences in the relative abundance of features between sample groupings. It is a compositional approach that makes no assumptions about feature distributions. However, it requires that all features have non-zero abundances so a pseudocount first needs to be added (1 is a typical pseudocount choice):

```         
qiime composition add-pseudocount \
  --i-table deblur_output/deblur_table_final.qza \
  --p-pseudocount 1 \
  --o-composition-table deblur_output/deblur_table_final_pseudocount.qza
```

Then ANCOM can be run with this command; note that CATEGORY is a placeholder for the text label of your category of interest from the metadata file.

```         
qiime composition ancom \
  --i-table deblur_output/deblur_table_final_pseudocount.qza \
  --m-metadata-file Blueberry_metadata_reduced.tsv \
  --m-metadata-column CATEGORY \
  --output-dir ancom_output
```

**Question 6: If you run more than one description column, you won't be able to overwrite the output directory. How can you put a second file in this directory? (Hint: you can find additional information at [QIIME 2 user documentation](https://docs.qiime2.org/2023.2/) or by typing** `qiime composition ancom --help`**)**

**Question 7: Does ANCOM identify any differentially abundant taxa between any of the metadata groups? If so, which one(s)?**

## *6. (Optional) Phylogenetic Robust Aitchison PCA*

New approaches to calculating beta diversity are still being invented. For example, Phylogenetic Robust Aitchison PCA (RPCA) aims to account for both the compositional and phyogenetic nature of microbiome data can be run using the `gemelli` toolbox. If you have time and want to explore this option, you can find instructions for installing `gemelli` [here](https://github.com/biocore/gemelli?tab=readme-ov-file) (you can install it directly into your QIIME2 environment), and follow the tutorial for running it [here](https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/Phylogenetic-RPCA-moving-pictures.ipynb). Installing and trying out new tools and methods is an important part of bioinformatics!

## Answers

**Question 1: What is a good rarefaction depth for diversity analysis?**

A cut-off of 4,000 reads will be sufficient: the curve plateaus around this depth and we won't exclude any samples.

**Question 2: are there any significant differences in alpha diversity between any of our metadata categories?**

There are significant differences in richness and phylogenetic diversity between forest and managed environment samples. There are no significant differences between bulk and rhizosphere soil.

**Question 3: which metadata category appears to provide more separation in the beta diversity PCoA plots?**

This is hard to say just by looking at the Emperor plots provided by QIIME2, but the forest/managed category appears to exhibit more distinct separation. The PERMANOVA test from the `beta-group-significance` command shows this as well.

**Question 4: what do you mean I'm getting an error message? It should be working fine, it worked last time! I can't see what I'm doing wro... oooooh there we go, I see now :)**

There is a typo in `--m-metadata-column Desription_1`: it's missing the "c" and should be `Description_1` ;)

**Question 5: can you identify any patterns between the metadata groups?**

Because stacked barcharts are limited in their analytical capabilities, it is hard to discern anything except very obvious patterns.

**Question 6: If you run more than one description column, you won't be able to overwrite the output directory. How can you put a second file in this directory? (Hint: you can find additional information at [QIIME 2 user documentation](https://docs.qiime2.org/2023.2/) or by typing** `qiime composition ancom --help`**)**

```         
qiime composition ancom \
  --i-table deblur_output/deblur_table_final_pseudocount.qza \
  --m-metadata-file Blueberry_metadata_reduced.tsv \
  --m-metadata-column Description_3 \
  --o-visualization ancom_output/description_4.qzv
```

**Question 7: Does ANCOM identify any differentially abundant taxa between any of the metadata groups? If so, which one(s)?**

ANCOM does not identify any taxa as differentially significant between bulk and rhizosphere soils, but three ASVs are identified as differentially abundant between forest and managed environments. You can look up the identity of each ASV in the taxonomy.tsv file.
