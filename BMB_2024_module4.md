---
layout: aws_tutorial_page
permalink: /BMB_2024_module4
title: AWS 2024 - Module 4
header1: Workshop Pages for Students
header2: Module 4 - Functional prediction and additional analyses
image: /site_images/BMB_2024_v1.png
home: https://bioinformaticsdotca.github.io/BMB_2024
---

# Module 4: Functional prediction and additional analyses

This tutorial is part of the 2024 Canadian Bioinformatics Workshops [Beginner Microbiome Analysis](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-workshop-2024-beginner) (St John's, NL, May 27-28).

**Author**: Robyn Wright

## Table of Contents

[Introduction](#introduction)\
[4.1 Prepare data from the end of Module 2](#41-prepare-data-from-the-end-of-module-2)\
[4.2 Start running PICRUSt2](#42-start-running-picrust2)\
[4.3 Read filtered output into R/Phyloseq](#43-read-filtered-output-into-rphyloseq)\
[4.4 Run MaAsLin2](#44-run-maaslin2)\
[4.5 Run ANCOM2](#45-run-ancom2)\
[4.6 Run ALDEx2](#46-run-aldex2)\
[4.7 Combine and plot differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_1](#47-combine-and-plot-differential-abundance-results-from-maaslin2-ancom2-and-aldex2-for-description_1)\
[4.8 Combine and plot differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_3](#48-combine-and-plot-differential-abundance-results-from-maaslin2-ancom2-and-aldex2-for-description_3)\
[4.9 Read PICRUSt2 output into R/Phyloseq](#49-read-picrust2-output-into-rphyloseq)\
[4.10 Look at PICRUSt2 alpha and beta diversity](#410-look-at-picrust2-alpha-and-beta-diversity)\
[4.11 Run PICRUSt2 differential abundance](#411-run-picrust2-differential-abundance)\
[4.12 Combine and plot PICRUSt2 differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_1](#412-combine-and-plot-picrust2-differential-abundance-results-from-maaslin2-ancom2-and-aldex2-for-description_1)\
[4.13 Combine and plot PICRUSt2 differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_3](#413-combine-and-plot-picrust2-differential-abundance-results-from-maaslin2-ancom2-and-aldex2-for-description_3)\
[Answers](#answers)

## Introduction

In this module, we'll be taking a look at two of the major ways in which we use amplicon sequencing data to learn about microbial communities; (1) predicting the functional capabilities of these communities based on the taxa that are present, and (2) identifying taxa that are significantly differentially abundant between different sample groups. 

There is some general information on functional prediction with PICRUSt2 and differential abundance testing with MaAsLin2, ANCOM2 and ALDEx2 below, but during the tutorial, you'll start by preparing and filtering the data that was obtained from QIIME2 at the end of Module 2, then will start PICRUSt2 running. As this may take a long time to run, you can run differential abundance testing on the filtered Module 2 outputs and examine these results before going back to look at the PICRUSt2 results. 

Throughout this module, there are some questions aimed to help your understanding of some of the key concepts. You'll find the [answers](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-4:-Functional-prediction-and-additional-analyses#answers) at the bottom of this page, but no one will be marking them. 

### Functional prediction with PICRUSt2 

PICRUSt2 is a tool that predicts the functional capacity of a microbial community based on the taxa that are present in amplicon sequencing data. The first iteration of [PICRUSt](https://www.nature.com/articles/nbt.2676) (**P**hylogenetic **I**nvestigation of **C**ommunities by **R**econstruction of **U**nobserved **St**ates) was developed by Morgan (and many other collaborators) during his postdoc, and this was expanded and improved upon to make [PICRUSt2](https://www.nature.com/articles/s41587-020-0548-6) by a previous PhD student in the lab, [Gavin Douglas](https://www.gavindouglas.ca/) (currently a postdoc at North Carolina State University). 

There are several other tools that have also been developed for this purpose, *e.g.*, [Tax4Fun2](https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/s40793-020-00358-7), [Piphillin](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6427-1) and [PanFP](https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-015-1462-8). I don't personally have experience including these, so we're going to be focussing on PICRUSt2, however, it's always important to understand the strengths and weaknesses of different bioinformatic tools before choosing which one to use for your own research. 

You can see full information on PICRUSt2 [here](https://github.com/picrust/picrust2/wiki), but it includes several key steps:
1. Alignment of ASVs to reference sequences using [HMMER](http://www.hmmer.org/)
2. Determining the optimal placement of ASVs into the reference tree using [EPA-NG](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy054/5079844) or and [GAPPA](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa070/5722201) to output a new tree incorporating the ASV placements (or [SEPP](https://www.worldscientific.com/doi/abs/10.1142/9789814366496_0024), to combine both of these tasks using less memory)
3. Inferring gene family copy numbers of ASVs using [castor](https://academic.oup.com/bioinformatics/article/34/6/1053/4582279)
4. Determining gene family copy numbers per sample
5. Inferring pathway abundances using [MinPath](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000465)

These steps can be run separately, but are all wrapped together into a single command for what we'll be running today. 

### Differential abundance testing with MaAsLin2, ANCOM2 and ALDEx2

Tens (if not more) of tools have been used for differential abundance testing of microbiome samples - identifying taxa that differ in abundance between sample groups (for example, treatment versus control) or with variables of interest (for example, correlations with pH or blood metabolite measurements). While these tools/methods have often been used interchangeably in the literature, there are large differences in the performance of different tools even on the same sample group. There have now been several large scale comparisons of different tools (*e.g.,* [Calgaro *et al.*](https://link.springer.com/article/10.1186/s13059-020-02104-1), [Thorsen *et al.*](https://link.springer.com/article/10.1186/s40168-016-0208-8)), and our lab also carried one out (published [here](https://www.nature.com/articles/s41467-022-28034-z)). While these studies have typically found that different methods lead to different biases in the results, as there is no one single *best* method for determining differentially abundant taxa, the best practice that we have for overcoming these issues currently is to use multiple differential abundance tests, report these clearly, and focus on the taxa that are identified by multiple tests. 

For these purposes, we have chose to use three tools, [MaAsLin2](https://huttenhower.sph.harvard.edu/maaslin/), [ANCOM2](https://www.nature.com/articles/s41467-020-17041-7) and [ALDEx2](https://www.biorxiv.org/content/10.1101/2023.10.21.563431v1). In our comparison of different tools, we found that both ANCOM2 and ALDEx2 controlled the false discovery rate well, but this came at the cost of reduced sensitivity, while MaAsLin2 had much better sensitivity than either ANCOM2 or ALDEx2, but this came at the cost of reduced specificity. In some of the recent studies from our lab, we have then chosen to focus on only the taxa that are identified by >two of these tools (*e.g.,* [this paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01662-3) or [this paper](https://www.nature.com/articles/s41598-024-60409-8)). 

In this tutorial, we'll be running each of these tools on the output from Module 2, and then we'll combine the results and look at the taxa that are identified. We'll then apply this to the PICRUSt2 output, too. 

## 4.1 Prepare data from the end of Module 2

First, we need to change directory to our workspace, create a new directory for this module, and change into that directory:
```
cd workspace
mkdir bmb_module4
cd bmb_module4/
```

Then, we'll create symlinks to the relevant outputs from Module 2, the denoised ASVs (deblur_output), the taxonomy classifications (taxa), and the metadata (Blueberry_metadata_reduced.tsv):
```
ln -s ~/workspace/bmb_module2/deblur_output/ .
ln -s ~/workspace/bmb_module2/taxa/ .
ln -s ~/workspace/bmb_module2/Blueberry_metadata_reduced.tsv .
```

Remember, we can use the ```ls``` command to check that everything copied across OK.

### Filter and export files

As PICRUSt2 can take quite a long time to run, we want to do some filtering on the files so that we'll have fewer ASVs. Usually we might investigate several different cutoffs for the minimum number of times an ASV should occur to be included, or the minimum prevalence of ASVs, but for the purposes of this workshop, we'll just use 50 and 2.

First, we need to reactivate the QIIME2 environment:
```
conda activate qiime2-amplicon-2024.2
```

Then we'll filter the ASVs based on the parameters I mentioned above:
```
qiime feature-table filter-features \
  --i-table deblur_output/deblur_table_final.qza \
  --p-min-frequency 50 \
  --p-min-samples 2 \
  --o-filtered-table final_table_filtered.qza
```

Now we'll filter the sequences to only include those that are in our new filtered table (```final_table_filtered.qza```):
```
qiime feature-table filter-seqs \
  --i-data deblur_output/representative_sequences.qza \
  --i-table final_table_filtered.qza \
  --o-filtered-data  representative_sequences_final_filtered.qza
```

Then we can export these files - first the feature table, which we'll export and then convert from .biom format to .txt format:
```
qiime tools export \
  --input-path final_table_filtered.qza \
  --output-path exports_filtered
  
biom convert \
  -i exports_filtered/feature-table.biom \
  -o exports_filtered/feature-table.txt \
  --to-tsv 
```

Then we can export the sequences file:
```
qiime tools export \
  --input-path representative_sequences_final_filtered.qza \
  --output-path exports_filtered
```
This will be exported as a .fasta file, which you can look at in the ```exports_filtered``` folder, if you like.

Finally, we'll export the taxa classifications that we have, we'll slightly modify the resulting file, and then deactivate the conda environment:
```
qiime tools export \
  --input-path taxa/classification.qza \
  --output-path taxa
  
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv
  
conda deactivate
```

## 4.2 Start running PICRUSt2

PICRUSt2 can take quite a long time to run - for PICRUSt2, as well as other programs that may take a while, there are several tools that are pre-installed on most Linux systems that we can use to make sure that our program carries on running even if we get disconnected from the server. One of the most frequently used ones is called ```tmux```. To activate it, just type in ```tmux``` and press enter. It should take a second to start up, and then load up with a similar looking command prompt to previously, but with a coloured bar at the bottom of the screen.

To get out of this window again, press ```ctrl```+```b``` at the same time, and then d. You should see your original command prompt and something like
```
[detached (from session 0)]
```

We can actually use tmux to have multiple sessions, so to see a list of the active sessions, use:
```
tmux ls
```

We can rename the tmux session that we just created with this:
```
tmux rename-session -t 0 picrust
```
Note that we know it was session 0 because it said that we detached from session 0 when we exited it.

If we want to re-enter this window, we use:
```tmux attach-session -t picrust```

Now, we can run PICRUSt2 inside this tmux session.

First, we'll activate the conda environment that has PICRUSt2 installed:
```
conda activate picrust2
```

Now, as we mentioned above, we'll just set PICRUSt2 to run for now, and we'll come back to see the output later:
```
picrust2_pipeline.py -s exports_filtered/dna-sequences.fasta -i exports_filtered/feature-table.biom -o picrust2_out_pipeline_filtered -p 4 -t sepp
```
You can see that here the options we're setting are:
- ```-s```: the fasta file of DNA sequences
- ```-i```: the feature table in .biom format
- ```-o```: the folder to save the PICRUSt2 output to
- ```-p```: the number of threads to use
- ```-t```: the method to use for placement of ASVs into the phylogenetic tree - note that the default is to use EPA-NG, but this takes much more memory, and as we are limited on these servers we are using SEPP

You can move onto the next sections of this module while PICRUSt2 runs. 

## 4.3 Read filtered output into R/Phyloseq

While PICRUSt2 runs, we're going to run the differential abundance section of this workshop.

Now, we'll go to RStudio on the server. Go to your browser and navigate to:
```
http://##.uhn-hpc.ca:8080
```
(where ## is the student number that you have been given).

The username is "ubuntu" and you'll be given the password in class. 

Everything that you're entering in here should be entered into the "Console" part of the page. This is like the Terminal that we've been using on the Server, but is uses the R programming language. Most statistical analyses in bioinformatics will use either R or Python - many bioinformaticians have a preference for one or the other, but most use both. In my opinion, each has different strengths and weaknesses - R has many built in programs and functions, while Python is much more customisable but perhaps has a steeper learning curve. The syntax varies slightly for both (and both are different from the Terminal that we have been using, which is the bash programming language), but there are also many similarities. RStudio is an Integrated Developing Environment (IDE) for R and Python - it is an application that provides further functionality than the R programming language itself. 

### Importing packages

First of all, we'll import the packages that we'll be using today:
```
library(reticulate)
use_condaenv('/home/ubuntu/CourseData/MIC_data/.conda/envs/r-env/')
library(phyloseq)
library(Maaslin2, lib.loc="/home/ubuntu/CourseData/MIC_data/.conda/envs/r-env/lib/R/library")
library(ANCOMBC, lib.loc="/home/ubuntu/CourseData/MIC_data/.conda/envs/r-env/lib/R/library")
library(microbiome, lib.loc="/home/ubuntu/CourseData/MIC_data/.conda/envs/r-env/lib/R/library")
library(ALDEx2)
library(ggplot2)
library(tidyr)
library(stringr)
library(vegan)
```
You'll notice that in each of these, we use the ```library()``` command, and inside the brackets we give the name of the package. In some cases, we also tell R where it can find the files associated with this package - this can be useful if you have multiple versions of something installed, or it is installed in a location that R doesn't know to search for by default. We've already installed these packages for you, but in most cases these could be installed by *e.g.,* ```install.packages("reticulate")```. In some cases, they might be installed using another package, BiocManager, *e.g.,* ```BiocManager::install("ALDEx2")```. If you're unsure, you can always google how to install a package and you can usually find the code to copy and paste, *e.g.,* "r install aldex2". You could also try typing ```??aldex2``` into the console and see what comes up. It should give you some information about that package.

### Setting variables

Next, we'll set the name of the folder where our results are saved:
```
folder = '/home/ubuntu/workspace/bmb_module4/'
```
Here, we have defined a variable called "folder" - we always need to be careful that we don't save over these variables by calling something else the same thing, but you can always check these in the "Environment" area in the top right where you should see that "folder" now exists. We could have called this anything, but it's always helpful to name our variables something that is simple and explains what they are, incase someone else is using our code. This is a string, as it is inside the punctuation marks ```' '``` or ```" "```, and this is used in programming to show that we are representing text rather than numbers. These variables can be really useful to avoid having to type in the same information multiple times, which always increases the chances that we'll make a mistake!

You can see how these would be different by running the following:
```{R}
random_string = "12345"
random_number = 12345
```
You should see how these are different by how they show in the "Environment" area in the top right. 

Now that we've explained the basics, we can read in some of our data.

### Read in the feature table

First, we'll read in the feature table. Note that we're calling it ```asv_table```:
```
asv_table <- read.csv(paste(folder, "exports_filtered/feature-table.txt", sep=""), sep='\t', skip=1, header=T)
```
Here we're combining a few different things. We're using the ```read.csv``` function - functions are modules of code that accomplish a specific task. We can write them ourselves, they are what is contained in the packages that we imported at the start, or there are several that are built in to R. If you want more information about what a function does as well as the input that it expects and the output that it gives, you can type in ```?read.csv```. 

We're also using the ```paste``` function to add together the full file path of the feature table. We're giving it the ```folder``` variable that we already defined, as well as the additional folder and file name (both as positional arguments), and finally the ```sep=""``` named argument - this tells the ```paste``` function that there should be no spaces or other punctuation between the two parts that are being pasted together. Try running just this part and see what the result is: ```paste(folder, "exports_filtered/feature-table.txt", sep="")``` - note that as we haven't given this a name, it is not saved as a variable. 

Next, we've given the ```sep='\t'``` variable - this time, we're telling the ```read.csv``` function that the file that we're importing is tab-delimited, then we've told it to skip the first line of the file ```skip=1```, and that there is a header ```header=T``` (where T is the same as TRUE). We've skipped the first line because this just reads "# Constructed from biom file" and is not actually data.

### Manipulate the feature table

We often find in programming that things are not quite in the format that is expected by the packages that we're using, so we'll make a few modifications to ```asv_table```:
```
asv_table_num = data.matrix(asv_table[,2:13]) #convert the ASV table to a numeric matrix
rownames(asv_table_num) = asv_table[,1] #give the matrix row names
```
You can see here that I've added ```#``` after the code and have written what each line does. You can use ```#``` to make comments through your code documents, to explain what each line is doing. Here you should also notice that we didn't save over the previous ```asv_table```, but created a new one called ```asv_table_num```. Now, we've just told R that it should be expecting numberic data within this table, and then we gave it the same row names as our previous ```asv_table``` object. 

### Read in the taxonomy information

Now, we'll be doing similar with our taxonomy information:
```
taxonomy <- read.csv(paste(folder, "taxa/taxonomy.tsv", sep=""), sep='\t')
```
If you print this out by typing ```taxonomy``` or by clicking on this in the "Environment" area, you'll see that we currently only have a single column that contains all of the taxonomy information. 

We want to split this so that we have columns for each of Domain, Phylum, Class, Order, Family, Genus and Species. Luckily, there is already a function that we can use for this:
```
taxonomy_split <- separate(data = taxonomy, col = taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\; ") #separate the taxonomy table so each phylogenetic level is its own column
```
You'll see that we tell the ```separate``` function that the data it should use is in the ```taxonomy``` table, the ```col``` (column) is called "taxonomy". The "into" named argument takes a list as input, which is defined in R with the ```c()```, and the ```taxonomy``` column should be split into a new column each time the ```;``` symbol is found. 

Again, you can take a look at the results in the "Environment" area. 

Now we'll do some final manipulations to get this into the right format:
```
taxonomy_split <- taxonomy_split[,-c(1,9)] #remove the Feature.ID and Confidence columns from the taxonomy table
rownames(taxonomy_split) <- taxonomy[,1] #and now give the taxonomy table the OTU IDs as row names
```

### Read in the metadata

Now we'll read in the metadata. I won't go through each stage step-by-step as hopefully you're getting the idea by now, but I've still added comments to each of the rows:
```
metadata <- read.csv(paste(folder, "Blueberry_metadata_reduced.tsv", sep=""), sep='\t')
samples = metadata[,2:3] #get the metadata columns
rownames(samples) = metadata[,1] #and add the sample names as row names
samples = data.frame(samples, stringsAsFactors = FALSE) #convert this to a data frame
```

### Combine these into a phyloseq object

[Phyloseq](https://joey711.github.io/phyloseq/) is a really useful package that contains many useful functions for analysing microbiome data. While it can be a bit fiddly to get our data into the format that it is expecting (most good packages should give you examples of how data should be formatted to go into them, although it does take a lot of practice to get good at quickly working out how this is different from what you have!), once that it is in this format, the analyses are then very easy to perform. 

```
ASV = otu_table(asv_table_num, taxa_are_rows = TRUE) #convert asv_table_num to an otu_table
TAX = tax_table(taxonomy_split) #convert taxonomy_split to a tax_table
taxa_names(TAX) <- rownames(taxonomy_split) #add names to TAX/the tax_table
SAMPLE = sample_data(samples) #convert samples to a sample_data
physeq = phyloseq(ASV, TAX, SAMPLE) #combine these all to make a phyloseq object
physeq #print this out to see what a phyloseq object looks like
```
Now, we have combined all of these different parts into one phyloseq object. Note that you can also import a phylogenetic tree to phyloseq, but we won't be using that for this part of the workshop. You can also get each of the individual objects back after performing manipulations, *e.g.,* ```otu_table(physeq)```.

As we typically find that ASVs are not always shared across many samples (and neither would we expect them to be, if they potentially represent species or strain level differences between taxa), we'll collapse the phyloseq object to the genus level. If you're not sure what the different taxonomy levels are called within the phyloseq object, you can always look at ```tax_table(physeq)```.

So we're collapsing at the genus level, which is "ta6":
```
physeq_genus = tax_glom(physeq, taxrank="ta6")
```
If you take a look at this taxonomy table (```tax_table(physeq_genus)```), you'll notice a lot of missing information (because many environmental ASVs don't have similar taxa within reference databases!), so we can add the full taxonomy information in to the ASV names, so that this is hopefully a little more informative:
```
all_tax = paste(tax_table(physeq_genus)[,2], tax_table(physeq_genus)[,3], tax_table(physeq_genus)[,4], tax_table(physeq_genus)[,5], tax_table(physeq_genus)[,6], sep=';')
taxa_names(physeq_genus) = all_tax
otu_table(physeq_genus)
```
Take a look at each of these if you want to see what the differences were in each step!

## 4.4 Run MaAsLin2

Now that we've prepared the phyloseq objects, it's time to run our first differential abundance test. Each of the differential abundance tools expects data to be normalised in a different way, and we're going to run MaAsLin2 with rarefied data.

So to do that, we'll need to rarefy it:
```
physeq_rare = rarefy_even_depth(physeq_genus, sample.size = min(sample_sums(physeq_genus)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
```
Remember that you can always use ```?rarefy_even_depth``` if you want to see what a function is doing!
Here we've rarefied all of the samples to the lowest number of reads, ASVs are replaced after sampling (```replace = TRUE```), those that are no longer present are trimmed (```trimOTUs = TRUE```) and the function will tell us about what it is doing (```verbose = TRUE```).

Next, we'll get the feature table and metadata that we're using:
```
feat_table = data.frame(t(otu_table(physeq_rare)), check.rows=T, check.names=T, stringsAsFactors=T)
metadata = data.frame(sample_data(physeq_rare), stringsAsFactors = F)
```
You'll notice that these will be formatted as tables again - we could have just read in some tables to be used with MaAsLin2, but it's useful to perform manipulations in Phyloseq first, and then to take these tables out again to ensure that all samples match up!

Now we'll run MaAsLin2, first with Description_1 and then with Description_3:
```
results_all <- Maaslin2(feat_table, metadata, paste(folder, "MaAsLin2_out_taxa_Description_1", sep=""), transform = "AST", fixed_effects = c("Description_1"), reference=paste("Description_1", "Bulk", sep=','), standardize = FALSE, plot_heatmap = T, plot_scatter = T)
results_all <- Maaslin2(feat_table, metadata, paste(folder, "MaAsLin2_out_taxa_Description_3", sep=""), transform = "AST", fixed_effects = c("Description_3"), reference=paste("Description_3", "Forest", sep=','), standardize = FALSE, plot_heatmap = T, plot_scatter = T)
```
Take a look at ```?Maaslin2``` to see what each of the inputs to this function are!

MaAsLin2 will save the results into a folder, so you can take a look in these folders and look at the ```all_results.tsv``` as well as ```significant_results.tsv``` files. You'll see that ```significant_results.tsv``` is just a subsection of ```all_results.tsv``` - have a look at one of these files. 

In the columns you'll see: (1) the genus ("feature"), (2) the metadata variable being tested ("metadata"), (3) for categorical features, the specific feature level for which the coefficient and significance of association is being reported ("value"), (4) the model effect size/coefficient value ("coef"), (5) the standard error from the model ("stderr"), (6) the number of samples used in the model for these association values ("N"), (7) the number of samples in which this feature is non-zero ("N.not.0"), (8) the nominal significance of this association ("pval"), and (9) the corrected significance  with p.adjust of this association ("qval"). 

You'll see that by default, MaAsLin2 reports any features with q-value <= 0.25 as being significant. 

## 4.5 Run ANCOM2

Now we'll run ANCOM2 - ANCOM2 performs normalisations within the package and it can also take a phyloseq object as input, so this is pretty straightforward to do, but it may take a minute or two to run:
```
ancom_out = ancombc(phyloseq=physeq_genus, formula="Description_1+Description_3", alpha=0.1)
```
You'll see that here we can add in both metadata variables at the same time, and if we added them as ```Description_1*Description_3``` then this would also look at the interaction between the variables, as well as each of the variables separately.

ANCOM does save some other variables to ```ancom_out```, but what tells us the most is in ```ancom_out$res``` - this gives us a load of different tables: ```ancom_out$res$beta```, ```ancom_out$res$se```, ```ancom_out$res$W```, ```ancom_out$res$p_val```, ```ancom_out$res$q_val``` and ```ancom_out$res$diff_abn```, where beta shows an indication of the effect size of the differences between groups, se is the standard errors of these differences, W is the test statistic, p_val is the p-values for tests and q_val is the adjusted p-values for tests. diff_abn indicates whether the tests were significant based on the alpha value given to the function (```alpha=0.1```). 

In this workshop, we'll just be making a new table that contains the W test statistics, p-values and q-values for tests, but if you're using this in your own research then you should think about whether you want to include some kind of effect size filter. 

Make the new table:
```
w = ancom_out$res$W
q = ancom_out$res$q_val
p = ancom_out$res$p_val
rownames(w) = w[,1] #rename the rows in w
w = w[,-c(1:2)] #remove columns 1-2 from w
rownames(q) = q[,1] #rename the rows in q
q = q[,-c(1:2)] #remove columns 1-2 from q
rownames(p) = p[,1] #rename the rows in p
p = p[,-c(1:2)] #remove columns 1-2 in p
colnames(w) = c("Description_1 W", "Description_3 W") #change the colnames in w
colnames(q) = c("Description_1 q", "Description_3 q") #change the colnames in q
colnames(p) = c("Description_1 p", "Description_3 p") #change the colnames in p
a_out = merge(w, q, by = 'row.names') #merge w and q based on their row names
rownames(a_out) = rownames(w) #rename the rows
a_out = a_out[,-1] #remove the duplicated column
a_out = merge(a_out, p, by = 'row.names') #merge a_out (containing w and q) with p based on their row names
rownames(a_out) = rownames(w) #rename the rows
a_out = a_out[,-1] #remove the duplicated column
write.csv(a_out, paste(folder, "ANCOM_taxa.csv", sep="")) #save the resulting table as a .csv file
```

## 4.6 Run ALDEx2

Finally, we'll run ALDEx2. Note that here the first step is to normalise the samples using a CLR normalisation, so you can see how all of the methods require similar steps - normalisation and then testing - but this is done in a slightly different way for each. We then use a Kruskal-Wallis test and a GLM ANOVA on the data, and we can do this test with multiple cores (```useMC=2```). You can see that we're again running this separately for each metadata variable, and then saving the outputs as .csv files:
```
x <- aldex.clr(otu_table(physeq_genus), sample_data(physeq_genus)$Description_1, mc.samples = 128, verbose=F, denom="all")
kw.test.d1 <- aldex.kw(x, useMC=2, verbose=FALSE)
write.csv(kw.test.d1, paste(folder, "ALDEx2_taxa_Description_1.csv", sep=""))

x <- aldex.clr(otu_table(physeq_genus), sample_data(physeq_genus)$Description_3, mc.samples = 128, verbose=F, denom="all")
kw.test.d3 <- aldex.kw(x, useMC=2, verbose=FALSE)
write.csv(kw.test.d3, paste(folder, "ALDEx2_taxa_Description_3.csv", sep=""))
```
If you take a look at these results files then you'll see that we have 4 columns: kw.ep, kw.eBH, glm.ep, glm.eBH - the p-values and Benjamini-Hochberg (BH)-corrected p-values for the Kruskal-Wallis and GLM ANOVA tests between the sample groups. There are several other options for testing within the ALDEx2 R package (including aldex.ttest and aldex.glm), but in this case, we can use the GLM ANOVA results for parametric tests and the Kruskal-Wallis results for non-parametric tests. We could also use aldex.effect to see the magnitude of the differences between groups. 

**Note**: parametric tests make assumptions about the distribution of the population that the samples are taken from (typically that it is normally distributed), while non-parametric tests are distribution free and can be used for non-normal variables. Microbiome data are not typically normally distributed, but we can test this. You don't need to run this next part, but can do if you like.

```{R}
library(stats)
shapiro.test(as.data.frame(otu_table(physeq_genus))$BB197)
physeq_genus_clr = microbiome::transform(physeq_genus, "clr")
shapiro.test(as.data.frame(otu_table(physeq_genus_clr))$BB197)
```
Transformations/normalisations of microbiome data are often used to make the usually non-normal microbiome data normal, however, in this case, it is not normal (significant p-value so the null hypothesis that the data are normal can be rejected). This means that we should use the Kruskal-Wallis test results. 

## 4.7 Combine and plot differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_1

In order to compare the results that we've got from the three differential abundance tools, we'll want to first import the results and do some manipulations on these tables. 

For MaAsLin2, we'll read in the table and then we'll rename some of the taxa names because some R formats don't like punctuation other than "_" or ".", so we'll convert these back to how they started (so they match up with the other tables):
```
maaslin = read.csv(paste(folder, "MaAsLin2_out_taxa_Description_1/significant_results.tsv", sep=""), sep="\t")
maaslin$feature = gsub("\\g__Burkholderia.Caballeronia.Paraburkholderia", "g__Burkholderia-Caballeronia-Paraburkholderia", maaslin$feature) #replace anything matching g__Burkholderia.Caballeronia.Paraburkholderia in the feature column with g__Burkholderia-Caballeronia-Paraburkholderia
maaslin$feature = gsub("\\RCP2.54", "RCP2-54", maaslin$feature) #replace anything matching RCP2.54 in the feature column with RCP2-54
maaslin$feature = gsub("\\.", ";", maaslin$feature) #replace any remaining . with ;
maaslin = maaslin[maaslin$qval <= 0.1, ]
```
If you now look at the ```maaslin``` object, you'll see that there's only one taxon that had a q-value below 0.1 remaining.

Now we'll do the same for ALDEx2 (although here we don't need to do the renaming):
```
aldex = read.csv(paste(folder, "ALDEx2_taxa_Description_1.csv", sep=""))
aldex = aldex[aldex$kw.eBH <= 0.1, ]
```
ALDEx2 didn't find any taxa to be significantly differentially abundant (BH-adjusted p-value <= 0.1), so this is an empty dataframe.

And then do the same again for ANCOM2:
```
ancom = read.csv(paste(folder, "ANCOM_taxa.csv", sep=""))
ancom = ancom[ancom$Description_1.q <= 0.1, ]
```
If we look at ```ancom```, then you'll see that there are actually quite a few taxa that ANCOM2 found to be significantly differentially abundant. 

I mentioned previously that often we will consider those taxa found by >=2 tests to be significantly differentially abundant as actually being differentially abundant. That isn't really possible to do in this case, because there is no overlap between the taxa identified by the different tests, but hopefully this demonstrates why it's important to look at the results from multiple tools. We'll carry on and have a look at these taxa anyway!

First we'll create a list of these taxa, and make sure that we have no duplicates (not really necessary here, but could be for other tests!):
```
taxa = c(maaslin$feature, aldex$X, ancom$X)
taxa = unique(taxa)
```

And now we'll convert the physeq_genus phyloseq object to relative abundance using the ```transform_sample_counts``` function:
```
physeq_relabun  = transform_sample_counts(physeq_genus, function(x) (x / sum(x))*100 )
```

### Plot DA single taxon

Next we can look at plotting the abundance of a single taxon. Have a look at the list of taxa (```taxa```). Replace the empty quotation marks below with one of these:
```
single_taxon = ""
```

And now we'll make a boxplot of it's abundance based on the categories in Description_1:
```
microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=single_taxon, line = NULL, violin = FALSE, na.rm = FALSE, show.points = TRUE) + ggtitle(str_replace(single_taxon, ';f__', '\nf__')) + ylab('Relative abundance (%)')
```
Note that if you didn't replace the quotation marks with the name of a taxon above, then you will get an error message! Otherwise, you should see the plot pop up in the "Plots" section on the bottom right of the window. 

To explain what we have just done here a little, we're using a function from the microbiome R package - sometimes there can be functions with the same name from multiple R packages, so adding ```microbiome::``` at the start tells R which one we are wanting to use. The first positional argument is the phyloseq object that we are using, and then everything else consists of named arguments: The metadata variable that we'd like to use for the x axis (and making the boxplots), what we'd like on the y axis (this needs to be a taxon that is in the phyloseq object), whether we'd like to map a variable onto the lines, whether we'd like this to be a violin plot rather than a boxplot, whether NAs should be removed, and whether points should be shown in addition to the boxes. 

This ```boxplot_abundance``` function builds upon the R package ggplot2, and we can therefore use additional arguments as we would with ggplot2. These are added with the + and are fairly self explanatory; we have added a title (ggtitle) and a y label (ylab). For the title, we've replaced part of the text with the "\n" which means that a line break has been added in to make this a little more readable. 

### Plot DA all taxa

We can also make the same plots within a for loop:
```
for(i in 1:length(taxa)) {
  print(microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=taxa[i]) + ggtitle(str_replace(taxa[i], ';f__', '\nf__')) + ylab('Relative abundance (%)'))
}
```
Now, if you press the forward and back buttons in the "Plots" pane, you can scroll through all of the taxa that we've identified as differentially abundant by one of the tests. 

#### For loops

For loops are an incredibly useful thing that we do in programming - to my knowledge, they exist in every programming language, and they allow us to repeat a section of code many times, with a different input each time. To explain what we are doing above, we can go through this step-by-step. If you already know about for loops, feel free to skip ahead to section 2.8. 

It is easiest to explain by showing you, so try running this:
```
print(taxa)
print(length(taxa))
```
You'll see that taxa is a list containing all of the taxa that we found to be significantly differentially abundant, and length(taxa) tells you how long it is, or how many taxa are in it.

Now try running this:
```
for(i in 1:length(taxa)) {
  print(i)
  print(taxa[i])
}
```
You should see that now, a number is printed out, and then that number item from the list is also printed out after it. The for loop is telling us that for every value of i between 1 and the length(taxa), print out i and then print out the i'th item in the list. You can use this for accessing particular items in lists, too, *e.g.*, ```taxa[2]```. 

This is a really simple example, but you can see how in the loop above, we just replaced ```single_taxon```, from when we just plotted one taxon, with ```taxa[i]``` within the loop, so that it would change each time. 

For loops can get incredibly complicated, can do lots of things, and can be nested inside other for loops, but we won't be doing anything more on them for now!

#### Back to the results

**Question**: Are there any taxa that you're surprised to see are significantly differentially abundant? Why?

## 4.8 Combine and plot differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_3

Now see if you can repeat these steps for Description_3! I've added the code from above here, but if you're struggling then the code you'll need is in the [answers](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-4:-Functional-prediction-and-additional-analyses#answers).

Import the MaAsLin2 results:
```
maaslin = read.csv(paste(folder, "MaAsLin2_out_taxa_Description_1/significant_results.tsv", sep=""), sep="\t")
maaslin$feature = gsub("\\g__Burkholderia.Caballeronia.Paraburkholderia", "g__Burkholderia-Caballeronia-Paraburkholderia", maaslin$feature) #replace anything matching g__Burkholderia.Caballeronia.Paraburkholderia in the feature column with g__Burkholderia-Caballeronia-Paraburkholderia
maaslin$feature = gsub("\\RCP2.54", "RCP2-54", maaslin$feature) #replace anything matching RCP2.54 in the feature column with RCP2-54
maaslin$feature = gsub("\\.", ";", maaslin$feature) #replace any remaining . with ;
maaslin = maaslin[maaslin$qval <= 0.1, ]
```
**Hint**: Make sure that you check the names that you have in the maaslin table! You might need to replace some other taxon names like we have for RCP2.54. 

Now we'll do the same for ALDEx2:
```
aldex = read.csv(paste(folder, "ALDEx2_taxa_Description_1.csv", sep=""))
aldex = aldex[aldex$kw.eBH <= 0.1, ]
```

And then do the same again for ANCOM2:
```
ancom = read.csv(paste(folder, "ANCOM_taxa.csv", sep=""))
ancom = ancom[ancom$Description_1.q <= 0.1, ]
```
**Hint**: This time it's not when reading in the file that you need to make changes!

Now recreate the list of taxa, but call it something different (*e.g.,* taxa_d3) so that it doesn't overwrite the existing list:
```
taxa = c(maaslin$feature, aldex$X, ancom$X)
taxa = unique(taxa)
```

We don't actually need to re-run this part, but it doesn't hurt if you do:
```
physeq_relabun  = transform_sample_counts(physeq_genus, function(x) (x / sum(x))*100 )
```

### Plot DA single taxon

Next we can look at plotting the abundance of a single taxon. Have a look at the list of taxa (```taxa```). Replace the empty quotation marks below with one of these:
```
single_taxon = ""
```

Now make the boxplot:
```
microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=single_taxon, line = NULL, violin = FALSE, na.rm = FALSE, show.points = TRUE) + ggtitle(str_replace(single_taxon, ';f__', '\nf__')) + ylab('Relative abundance (%)')
```
**Hint**: Just replace the "Description_1" with "Description_3" here. 

### Plot DA all taxa

Now make all of the boxplots:
```
for(i in 1:length(taxa)) {
  print(microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=taxa[i]) + ggtitle(str_replace(taxa[i], ';f__', '\nf__')) + ylab('Relative abundance (%)'))
}
```

## 4.9 Read PICRUSt2 output into R/Phyloseq

PICRUSt2 should have finished running while you were looking at the differential abundance results, so now we can have a look at those results. To check everything, go back to the Terminal/Putty window and re-attach to the tmux session that you were using:
```
tmux attach-session -t picrust
```
You should hopefully see that it is completed. 

There are a few files that we'll take a look at now, so to prepare for that, we'll just unzip them:
```
gunzip picrust2_out_pipeline_filtered/*.gz
```

There is a key file that we can take a look at to see how well PICRUSt2 is likely to work for our data, and that is the ```marker_predicted_and_nsti.tsv``` file. Download this and open it with Excel (or similar). You'll see that it has three columns - the first contains the ASV name, the second the number of 16S rRNA gene copies that this ASV is predicted to have, and the third the NSTI. The NSTI is the Nearest Sequenced Taxon Index and refers to the distance within the phylogenetic tree of that ASV to the closest relative that it has in the reference database. A value of 0 would indicate that the database that PICRUSt2 uses contains a genome with an identical sequence. By default, PICRUSt2 excludes all ASVs with a value of above 2, but the lower these NSTI values are, the better the predictions are likely to be. In our case, the median NSTI is 0.171. This is not 0, but it is at least a long way off 2.

In this folder, you'll also see ```EC_predicted.tsv``` and ```KO_predicted.tsv``` - both of these show the number of copies of Enzyme Commission (EC) numbers and KEGG orthologs (KO), respectively, predicted to be within each ASV. You'll also see ```out.tre``` which contains a tree of the reference as well as study 16S sequences. 

The ```intermediate``` folder shows us some of the intermediate files that were produced/used, but we don't really need those for now.

Finally, the files that we are likely most interested in are ```picrust2_out_pipeline_filtered/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz```, ```picrust2_out_pipeline_filtered/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz``` and ```picrust2_out_pipeline_filtered/pathways_out/path_abun_unstrat.tsv.gz```.

We're going to focus on the pathways file today, because it groups the other functions into functional categories, but this could be repeated in the same way with either the KO or EC results. 

Note that if we were interested in seeing which ASVs contribute to the functions within each sample then we would include the ```--stratified``` option when running PICRUSt2, although this increases the time taken to run. 

So now try to do the same for the PICRUSt2 results as you did for the taxa results above in RStudio. I've again copied that code below here to avoid you needing to scroll up and down, and the code needed is in the [answers](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-4:-Functional-prediction-and-additional-analyses#answers) if you are struggling.

Read in the feature table:
```
asv_table <- read.csv(paste(folder, "exports_filtered/feature-table.txt", sep=""), sep='\t', skip=1, header=T)
```
**Hint**: You'll need to unzip the ```picrust2_out_pipeline_filtered/pathways_out/path_abun_unstrat.tsv.gz``` file before you can start.
**Hint**: You'll also need to change the file path ahead of the file name, and you won't want to skip the first line of the file anymore!
You will probably also want to name this something 

Manipulate the feature table:
```
asv_table_num = data.matrix(asv_table[,2:13]) #convert the ASV table to a numeric matrix
rownames(asv_table_num) = asv_table[,1] #give the matrix row names
```
**Hint**: Remember to change the table/object names everywhere that they occur!

We don't have taxonomy information this time, so that part isn't necessary, but next we'll want to read in the metadata. This will actually be exactly the same as before, but it won't hurt if you do it again.
```
metadata <- read.csv(paste(folder, "Blueberry_metadata_reduced.tsv", sep=""), sep='\t')
samples = metadata[,2:3] #get the metadata columns
rownames(samples) = metadata[,1] #and add the sample names as row names
samples = data.frame(samples, stringsAsFactors = FALSE) #convert this to a data frame
```

Now combine everything into the phyloseq object:
```
ASV = otu_table(asv_table_num, taxa_are_rows = TRUE) #convert asv_table_num to an otu_table
TAX = tax_table(taxonomy_split) #convert taxonomy_split to a tax_table
taxa_names(TAX) <- rownames(taxonomy_split) #add names to TAX/the tax_table
SAMPLE = sample_data(samples) #convert samples to a sample_data
physeq = phyloseq(ASV, TAX, SAMPLE) #combine these all to make a phyloseq object
physeq #print this out to see what a phyloseq object looks like
```
**Hint**: Make sure that you remember to change the names of everything and that you won't need the taxonomy information!
We also aren't collapsing this at the genus level because there is no taxonomy information/higher levels to collapse on.

Now, we can convert the pathways to relative abundance, take the top 100 pathways for differential abundance testing, and convert the numbers into integers (whole numbers) for some of the alpha diversity calculations. 
```
physeq_pwy_relabun = transform_sample_counts(physeq_pwy, function(x) (x / sum(x))*100 )
top_100_abun <- names(sort(taxa_sums(physeq_pwy_relabun), TRUE)[1:100]) #get most abundant pathways
physeq_pwy_top_100 = prune_taxa(top_100_abun, physeq_pwy) #now filter the table to have only the most abundant pathways
otu_table(physeq_pwy_top_100) = round(otu_table(physeq_pwy_top_100)) #round the abundance table so that everything is whole numbers
mode(pwy_table_num) = "integer" #convert the table to integers (whole numbers)
PWY_int = otu_table(pwy_table_num, taxa_are_rows = TRUE)
physeq_pwy_int = phyloseq(PWY_int, SAMPLE) #make a new phyloseq object with this
```
**Note**: These above commands assume that you called your pathway phyloseq object physeq_pwy and your pathway equivalent of asv_table_num pwy_table_num. You'll need to change this accordingly for whatever you called these. 

## 4.10 Look at PICRUSt2 alpha and beta diversity

Now we will just take a quick look at these to see how they compare with the analyses that we did in Module 3 based on the taxonomy data. 

First we will plot the Alpha diversity:
```{r}
plot_richness(physeq_pwy_int, x="Description_1", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
plot_richness(physeq_pwy_int, x="Description_3", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
```
Luckily, phyloseq has some built-in functions for these that make it very easy. You can see more information about these functions with ```?plot_richness```.

Then have a look at the Beta diversity:
```{r}
ps.ord <- ordinate(physeq_pwy_relabun, "PCoA", "bray")
plot_ordination(physeq_pwy_relabun, ps.ord, type="samples", color="Description_1", shape="Description_3") 
```
You should notice that we've chosen to colour the samples by Description_1 and the shape is based on Description_3. 

**Question**: Can you see differences between the samples based on the metadata variables? 
**Question**: What are these differences like compared with those seen in the taxonomy data? 

## 4.11 Run PICRUSt2 differential abundance

Now we'll try running the differential abundance tests. 

First, I'll help out with rarefying the pathway data:
```
physeq_rare_pwy = rarefy_even_depth(physeq_pwy, sample.size = min(sample_sums(physeq_pwy)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
physeq_rare_pwy_top_100 = prune_taxa(top_100_abun, physeq_rare_pwy)
```

But then I've again, I've copied in the code that we used above but you should modify it to work with your pathway phyloseq objects. Code that works is in the [answers](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-4:-Functional-prediction-and-additional-analyses#answers). 

Get the feature tables and metadata:
```
feat_table = data.frame(t(otu_table(physeq_rare)), check.rows=T, check.names=T, stringsAsFactors=T)
metadata = data.frame(sample_data(physeq_rare), stringsAsFactors = F)
```

Now we'll run MaAsLin2, first with Description_1 and then with Description_3:
```
results_all <- Maaslin2(feat_table, metadata, paste(folder, "MaAsLin2_out_taxa_Description_1", sep=""), transform = "AST", fixed_effects = c("Description_1"), reference=paste("Description_1", "Bulk", sep=','), standardize = FALSE, plot_heatmap = T, plot_scatter = T)
results_all <- Maaslin2(feat_table, metadata, paste(folder, "MaAsLin2_out_taxa_Description_3", sep=""), transform = "AST", fixed_effects = c("Description_3"), reference=paste("Description_3", "Forest", sep=','), standardize = FALSE, plot_heatmap = T, plot_scatter = T)
```
**Hint**: Make sure that you remember to change the names of the folders (*e.g.,* switch "taxa" for "pathway" to make sure that you don't write over the previous results!)

Run ANCOM2:
```
ancom_out = ancombc(phyloseq=physeq_genus, formula="Description_1+Description_3", alpha=0.1)
```

Make the table of ANCOM results:
```
w = ancom_out$res$W
q = ancom_out$res$q_val
p = ancom_out$res$p_val
rownames(w) = w[,1] #rename the rows in w
w = w[,-c(1:2)] #remove columns 1-2 from w
rownames(q) = q[,1] #rename the rows in q
q = q[,-c(1:2)] #remove columns 1-2 from q
rownames(p) = p[,1] #rename the rows in p
p = p[,-c(1:2)] #remove columns 1-2 in p
colnames(w) = c("Description_1 W", "Description_3 W") #change the colnames in w
colnames(q) = c("Description_1 q", "Description_3 q") #change the colnames in q
colnames(p) = c("Description_1 p", "Description_3 p") #change the colnames in p
a_out = merge(w, q, by = 'row.names') #merge w and q based on their row names
rownames(a_out) = rownames(w) #rename the rows
a_out = a_out[,-1] #remove the duplicated column
a_out = merge(a_out, p, by = 'row.names') #merge a_out (containing w and q) with p based on their row names
rownames(a_out) = rownames(w) #rename the rows
a_out = a_out[,-1] #remove the duplicated column
write.csv(a_out, paste(folder, "ANCOM_taxa.csv", sep="")) #save the resulting table as a .csv file
```
**Hint**: It's not really necessary to change everything here, just make sure that you change the name of the final output file!

Run ALDEx2:
```
x <- aldex.clr(otu_table(physeq_genus), sample_data(physeq_genus)$Description_1, mc.samples = 128, verbose=F, denom="all")
kw.test.d1 <- aldex.kw(x, useMC=2, verbose=FALSE)
write.csv(kw.test.d1, paste(folder, "ALDEx2_taxa_Description_1.csv", sep=""))

x <- aldex.clr(otu_table(physeq_genus), sample_data(physeq_genus)$Description_3, mc.samples = 128, verbose=F, denom="all")
kw.test.d3 <- aldex.kw(x, useMC=2, verbose=FALSE)
write.csv(kw.test.d3, paste(folder, "ALDEx2_taxa_Description_3.csv", sep=""))
```
**Hint**: Make sure again that you change the names of the output files and not just the input phyloseq objects!

## 4.12 Combine and plot PICRUSt2 differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_1

Again, I'll give the code that was used above, but you'll need to modify it. Code that works can be found in the [answers](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-4:-Functional-prediction-and-additional-analyses#answers).

Import the MaAsLin2 results:
```
maaslin = read.csv(paste(folder, "MaAsLin2_out_taxa_Description_1/significant_results.tsv", sep=""), sep="\t")
maaslin$feature = gsub("\\g__Burkholderia.Caballeronia.Paraburkholderia", "g__Burkholderia-Caballeronia-Paraburkholderia", maaslin$feature) #replace anything matching g__Burkholderia.Caballeronia.Paraburkholderia in the feature column with g__Burkholderia-Caballeronia-Paraburkholderia
maaslin$feature = gsub("\\RCP2.54", "RCP2-54", maaslin$feature) #replace anything matching RCP2.54 in the feature column with RCP2-54
maaslin$feature = gsub("\\.", ";", maaslin$feature) #replace any remaining . with ;
maaslin = maaslin[maaslin$qval <= 0.1, ]
```
**Hint**: It should just be within the ```read.csv``` function that you need to make changes! You probably won't need the changes to maaslin$feature that were being made previously, but you should check whethere there are any other changes that need to be made, and it won't hurt to have this left in there. 

Now we'll do the same for ALDEx2:
```
aldex = read.csv(paste(folder, "ALDEx2_taxa_Description_1.csv", sep=""))
aldex = aldex[aldex$kw.eBH <= 0.1, ]
```

And then do the same again for ANCOM2:
```
ancom = read.csv(paste(folder, "ANCOM_taxa.csv", sep=""))
ancom = ancom[ancom$Description_1.q <= 0.1, ]
```

Now make the list of pathways, but call it something different (*e.g.,* pwys) so that it doesn't overwrite the existing list:
```
taxa = c(maaslin$feature, aldex$X, ancom$X)
taxa = unique(taxa)
```

### Plot DA single pathway

Next we can look at plotting the abundance of a single pathways. Have a look at the list of pathways (```pwys```). Replace the empty quotation marks below with one of these:
```
single_taxon = ""
```
**Note**: It's probably good to change this to something like "single_pwy", although it's not necessary for actually making this work - just so that this is clear for the future!

Now make the boxplot:
```
microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=single_taxon, line = NULL, violin = FALSE, na.rm = FALSE, show.points = TRUE) + ggtitle(str_replace(single_taxon, ';f__', '\nf__')) + ylab('Relative abundance (%)')
```

### Plot DA all pathways

Now make all of the boxplots:
```
for(i in 1:length(taxa)) {
  print(microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=taxa[i]) + ggtitle(str_replace(taxa[i], ';f__', '\nf__')) + ylab('Relative abundance (%)'))
}
```

## 4.13 Combine and plot PICRUSt2 differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_3

Now do the same for Description_3! I think that you know the drill by now...

Import the MaAsLin2 results:
```
maaslin = read.csv(paste(folder, "MaAsLin2_out_taxa_Description_1/significant_results.tsv", sep=""), sep="\t")
maaslin$feature = gsub("\\g__Burkholderia.Caballeronia.Paraburkholderia", "g__Burkholderia-Caballeronia-Paraburkholderia", maaslin$feature) #replace anything matching g__Burkholderia.Caballeronia.Paraburkholderia in the feature column with g__Burkholderia-Caballeronia-Paraburkholderia
maaslin$feature = gsub("\\RCP2.54", "RCP2-54", maaslin$feature) #replace anything matching RCP2.54 in the feature column with RCP2-54
maaslin$feature = gsub("\\.", ";", maaslin$feature) #replace any remaining . with ;
maaslin = maaslin[maaslin$qval <= 0.1, ]
```
**Hint**: It should just be within the ```read.csv``` function that you need to make changes! You probably won't need the changes to maaslin$feature that were being made previously, but you should check whethere there are any other changes that need to be made, and it won't hurt to have this left in there. 

Now we'll do the same for ALDEx2:
```
aldex = read.csv(paste(folder, "ALDEx2_taxa_Description_1.csv", sep=""))
aldex = aldex[aldex$kw.eBH <= 0.1, ]
```

And then do the same again for ANCOM2:
```
ancom = read.csv(paste(folder, "ANCOM_taxa.csv", sep=""))
ancom = ancom[ancom$Description_1.q <= 0.1, ]
```

Now make the list of pathways, but call it something different (*e.g.,* pwys) so that it doesn't overwrite the existing list:
```
taxa = c(maaslin$feature, aldex$X, ancom$X)
taxa = unique(taxa)
```

### Plot DA single pathway

Next we can look at plotting the abundance of a single pathways. Have a look at the list of pathways (```pwys```). Replace the empty quotation marks below with one of these:
```
single_taxon = ""
```
**Note**: It's probably good to change this to something like "single_pwy", although it's not necessary for actually making this work - just so that this is clear for the future!

Now make the boxplot:
```
microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=single_taxon, line = NULL, violin = FALSE, na.rm = FALSE, show.points = TRUE) + ggtitle(str_replace(single_taxon, ';f__', '\nf__')) + ylab('Relative abundance (%)')
```

### Plot DA all pathways

Now make all of the boxplots:
```
for(i in 1:length(taxa)) {
  print(microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=taxa[i]) + ggtitle(str_replace(taxa[i], ';f__', '\nf__')) + ylab('Relative abundance (%)'))
}
```

# Answers

## Answers 4.8 Combine and plot differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_3

Read in MaAsLin2 results:
```
maaslin = read.csv(paste(folder, "MaAsLin2_out_taxa_Description_3/significant_results.tsv", sep=""), sep="\t")
maaslin$feature = gsub("\\g__Burkholderia.Caballeronia.Paraburkholderia", "g__Burkholderia-Caballeronia-Paraburkholderia", maaslin$feature)
maaslin$feature = gsub("\\RCP2.54", "RCP2-54", maaslin$feature)
maaslin$feature = gsub("\\__WPS.2", "__WPS-2", maaslin$feature)
maaslin$feature = gsub("\\JG36.TzT.191", "JG36-TzT-191", maaslin$feature)
maaslin$feature = gsub("\\.", ";", maaslin$feature)
maaslin = maaslin[maaslin$qval <= 0.1, ]
```

Read in ALDEx2 results:
```
aldex = read.csv(paste(folder, "ALDEx2_taxa_Description_3.csv", sep=""))
aldex = aldex[aldex$kw.eBH <= 0.1, ]
```

Read in ANCOM2 results:
```
ancom = read.csv(paste(folder, "ANCOM_taxa.csv", sep=""))
ancom = ancom[ancom$Description_3.q <= 0.1, ]
```

Combine the results:
```
taxa_d3 = c(maaslin$feature, aldex$X, ancom$X)
taxa_d3 = unique(taxa_d3)
```

Plot a single taxon:
```
single_taxon = ""
#e.g. single_taxon="p__Actinobacteriota;c__Actinobacteria;o__Catenulisporales;f__Actinospicaceae;g__Actinospica"

microbiome::boxplot_abundance(physeq_relabun, x='Description_3', y=single_taxon, line = NULL, violin = FALSE, na.rm = FALSE, show.points = TRUE) + ggtitle(str_replace(single_taxon, ';f__', '\nf__')) + ylab('Relative abundance (%)')
```

Plot all taxa:
```
for (i in 1:length(taxa_d3)) {
  print(microbiome::boxplot_abundance(physeq_relabun, x='Description_1', y=taxa_d3[i]) + ggtitle(str_replace(taxa_d3[i], ';f__', '\nf__')) + ylab('Relative abundance (%)'))
}
```

## Answers 4.9 Read PICRUSt2 output into R/Phyloseq

Unzip the pathways file:
```
gunzip picrust2_out_pipeline_filtered/pathways_out/path_abun_unstrat.tsv.gz
```

Read in the pathways table to R:
```
folder = '/home/ubuntu/workspace/bmb_module4/'

pwy_table <- read.csv(paste(folder, "picrust2_out_pipeline_filtered/pathways_out/path_abun_unstrat.tsv", sep=""), sep='\t', header=T)
pwy_table_num = data.matrix(pwy_table[,2:13]) #convert the ASV table to a numeric matrix
rownames(pwy_table_num) = pwy_table[,1] #give the matrix row names
```

Read in the metadata to R:
```
metadata <- read.csv(paste(folder, "Blueberry_metadata_reduced.tsv", sep=""), sep='\t')
samples = metadata[,2:3] #get the metadata columns
rownames(samples) = metadata[,1] #and add the sample names as row names
samples = data.frame(samples, stringsAsFactors = FALSE) #convert this to a data frame
```

Make the phyloseq object:
```
PWY = otu_table(pwy_table_num, taxa_are_rows = TRUE)
SAMPLE = sample_data(samples)
physeq_pwy = phyloseq(PWY, SAMPLE)
```

Convert the pathways to relative abundance, get the top 100, and convert the file to integers:
```
physeq_pwy_relabun = transform_sample_counts(physeq_pwy, function(x) (x / sum(x))*100 )
top_100_abun <- names(sort(taxa_sums(physeq_pwy_relabun), TRUE)[1:100]) #get most abundant pathways
physeq_pwy_top_100 = prune_taxa(top_100_abun, physeq_pwy) #now filter the table to have only the most abundant pathways
otu_table(physeq_pwy_top_100) = round(otu_table(physeq_pwy_top_100)) #round the abundance table so that everything is whole numbers
mode(pwy_table_num) = "integer" #convert the table to integers (whole numbers)
PWY_int = otu_table(pwy_table_num, taxa_are_rows = TRUE)
physeq_pwy_int = phyloseq(PWY_int, SAMPLE) #make a new phyloseq object with this
```
**Note**: This is also shown above, but I've put it all here again so that it's easier to follow.

## Answers 4.11 Run PICRUSt2 differential abundance

Run MaAsLIn2:
```
feat_table = data.frame(t(otu_table(physeq_pwy_top_100)), check.rows=F, check.names=F, stringsAsFactors=F)
metadata = data.frame(sample_data(physeq_pwy_top_100), stringsAsFactors = F)
results_all <- Maaslin2(feat_table, metadata, paste(folder, "MaAsLin2_out_pathways_Description_1", sep=""), transform = "AST", fixed_effects = c("Description_1"), reference=paste("Description_1", "Bulk", sep=','), standardize = FALSE, plot_heatmap = T, plot_scatter = T)
results_all <- Maaslin2(feat_table, metadata, paste(folder, "MaAsLin2_out_pathways_Description_3", sep=""), transform = "AST", fixed_effects = c("Description_3"), reference=paste("Description_3", "Forest", sep=','), standardize = FALSE, plot_heatmap = T, plot_scatter = T)
```

Run ALDEx2:
```
x <- aldex.clr(otu_table(physeq_pwy_top_100), sample_data(physeq_pwy_top_100)$Description_1, mc.samples = 128, verbose=F, denom="all")
kw.test.d1 <- aldex.kw(x, useMC=2, verbose=FALSE)
write.csv(kw.test.d1, paste(folder, "ALDEx2_pathways_Description_1.csv", sep=""))

x <- aldex.clr(otu_table(physeq_pwy_top_100), sample_data(physeq_pwy_top_100)$Description_3, mc.samples = 128, verbose=F, denom="all")
kw.test.d3 <- aldex.kw(x, useMC=2, verbose=FALSE)
write.csv(kw.test.d3, paste(folder, "ALDEx2_pathways_Description_3.csv", sep=""))
```

ANCOM2:
```
ancom_out = ancombc(phyloseq=physeq_pwy_top_100, formula="Description_1+Description_3", alpha=0.1)
w = ancom_out$res$W
q = ancom_out$res$q_val
p = ancom_out$res$p_val
rownames(w) = w[,1]
w = w[,-c(1:2)]
rownames(q) = q[,1]
q = q[,-c(1:2)]
rownames(p) = p[,1]
p = p[,-c(1:2)]
colnames(w) = c("Description_1 W", "Description_3 W")
colnames(q) = c("Description_1 q", "Description_3 q")
colnames(p) = c("Description_1 p", "Description_3 p")
a_out = merge(w, q, by = 'row.names')
rownames(a_out) = rownames(w)
a_out = a_out[,-1]
a_out = merge(a_out, p, by = 'row.names')
rownames(a_out) = rownames(w)
a_out = a_out[,-1]
write.csv(a_out, paste(folder, "ANCOM_pathways.csv", sep=""))
```

## Answers 4.12 Combine and plot PICRUSt2 differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_1

Read in MaAsLin2 results:
```
maaslin = read.csv(paste(folder, "MaAsLin2_out_pathways_Description_1/significant_results.tsv", sep=""), sep="\t")
maaslin$feature = gsub("\\.", "-", maaslin$feature)
maaslin = maaslin[maaslin$qval <= 0.1, ]
```

Read in ALDEx2 results:
```
aldex = read.csv(paste(folder, "ALDEx2_pathways_Description_1.csv", sep=""))
aldex = aldex[aldex$kw.eBH <= 0.1, ]
```

Read in ANCOM2 results:
```
ancom = read.csv(paste(folder, "ANCOM_pathways.csv", sep=""))
ancom = ancom[ancom$Description_1.q <= 0.1, ]
```

Combine results together and get the unique pathways:
```
pathways = c(ancom$X, aldex$X, maaslin$feature)
pathways = unique(pathways)
```

Plot a single pathway:
```
single_pathway = ""
#e.g. single_pathway="PWY-7222"

microbiome::boxplot_abundance(physeq_pwy_relabun, x='Description_1', y=single_pathway, line = NULL, violin = FALSE, na.rm = FALSE, show.points = TRUE) + ggtitle(single_pathway) + ylab('Relative abundance (%)')
```

Plot all pathways:
```
for(i in 1:length(pathways)) {
  print(microbiome::boxplot_abundance(physeq_pwy_relabun, x='Description_1', y=pathways[i]) + ggtitle(pathways[i]) + ylab('Relative abundance (%)'))
}
```

## Answers 4.13 Combine and plot PICRUSt2 differential abundance results from MaAsLin2, ANCOM2 and ALDEx2 for Description_3

Read in MaAsLin2 results:
```
maaslin = read.csv(paste(folder, "MaAsLin2_out_pathways_Description_3/significant_results.tsv", sep=""), sep="\t")
maaslin$feature = gsub("\\.", "-", maaslin$feature)
maaslin = maaslin[maaslin$qval <= 0.1, ]
```

Read in ALDEx2 results:
```
aldex = read.csv(paste(folder, "ALDEx2_pathways_Description_3.csv", sep=""))
aldex = aldex[aldex$kw.eBH <= 0.1, ]
```

Read in ANCOM2 results:
```
ancom = read.csv(paste(folder, "ANCOM_pathways.csv", sep=""))
ancom = ancom[ancom$Description_3.q <= 0.1, ]
```

Combine results together and get the unique pathways:
```
pathways_d3 = c(ancom$X, aldex$X, maaslin$feature)
pathways_d3 = unique(pathways_d3)
```

Plot a single pathway:
```
single_pathway = ""
#e.g. single_pathway="PWY-7222"

microbiome::boxplot_abundance(physeq_pwy_relabun, x='Description_3', y=single_pathway, line = NULL, violin = FALSE, na.rm = FALSE, show.points = TRUE) + ggtitle(single_pathway) + ylab('Relative abundance (%)')
```

Plot all pathways:
```
for(i in 1:length(pathways_d3)) {
  print(microbiome::boxplot_abundance(physeq_pwy_relabun, x='Description_3', y=pathways_d3[i]) + ggtitle(pathways_d3[i]) + ylab('Relative abundance (%)'))
}
```
