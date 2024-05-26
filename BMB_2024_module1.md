---
layout: aws_tutorial_page
permalink: /BMB_2024_module1
title: AWS 2024 - Module 1
header1: Workshop Pages for Students
header2: Module 1 - Introduction to sequencing data analysis
image: /site_images/BMB_2024_v1.png
home: https://bioinformaticsdotca.github.io/BMB_2024
---

# Module 1: Introduction to sequencing data analysis

This tutorial is part of the 2024 Canadian Bioinformatics Workshops [Beginner Microbiome Analysis](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-workshop-2024-beginner) (St John's, NL, May 27-28).

**Author**: Robyn Wright

## Table of Contents

[Introduction](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#introduction)\
[1.1 Log into your AWS instance](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#11-log-into-your-aws-instance)\
[1.2 Creating directories and moving around on the command line](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#12-creating-directories-and-moving-around-on-the-command-line)\
[1.3 Use wget to download files](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#13-use-wget-to-download-files)\
[1.4 Move these files to the directories](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#14-move-these-files-to-the-directories)\
[1.5 Zip and unzip these files](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#15-zip-and-unzip-these-files)\
[1.6 Create new tar archive of files](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#16-create-new-tar-archive-of-files)\
[1.7 Unzip tar archive](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#17-unzip-tar-archive)\
[1.8 Look at fasta and fastq files with less](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#18-look-at-fasta-and-fastq-files-with-less)\
[1.9 Installing programs to the server](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#19-installing-programs-to-the-server)\
[1.10 Conda environments](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#110-conda-environments)\
[1.11 Install fastqc and multiqc](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#111-install-fastqc-and-multiqc)\
[1.12 Perform quality control on fastq files](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#112-perform-quality-control-on-fastq-files)\
[1.13 htop - looking at the number of processes we have available or running](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#113-htop---looking-at-the-number-of-processes-we-have-available-or-running)\
[1.14 Back to the quality control](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#114-back-to-the-quality-control)\
[1.15 Quality control on the data that we'll be using](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#115-quality-control-on-the-data-that-well-be-using)\
[1.16 Install QIIME2](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#116-install-qiime2)\
[Answers](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024-Beginner-Module-1:-Introduction-to-sequencing-data-analysis#answers)

## Introduction

In this module, we're going to be learning the basics of moving around on the command line, creating directories, zipping and unzipping folders, making tar archives and installing programs on a server. We'll also be downloading files to the server, investigating the format of fasta and fastq files, and looking at the quality of sequenced samples. 

Throughout this module, there are some questions aimed to help your understanding of some of the key concepts. You'll find the answers at the bottom of this page, but no one will be marking them. 

## 1.1 Log into your AWS instance

Hopefully you were able to attend the session on Friday where we went over this, but I've added the instructions below as a reminder. You can also follow the instructions [here](https://github.com/bioinformaticsdotca/AWS_stuff/blob/master/Logging_into_the_Amazon_cloud.md).

### Mac/Linux

First off, you'll need to save the .pem key file that you've been given. 

Next, you'll need to change directory to where your key file is stored. Locate the file in your Finder, right click the file, hold the option key and click on "Copy BMB.pem as pathname". 

Then, go to your Terminal window (search for "Terminal" in Launchpad).

You can just remove the BMB.pem from the end to get the file path. E.g., my full path name that I copied is ```/Users/robynwright/Dropbox/Langille_Lab_postdoc/CBW_2024/student_test/BMB.pem```, so:
```
cd /Users/robynwright/Dropbox/Langille_Lab_postdoc/CBW_2024/student_test/
```

Then run:
```
chmod 600 BMB.pem
```
Make sure that you replace BMB.pem with whatever the file name is that you saved it as!

Now we can connect to the server:
```
ssh -i BMB.pem ubuntu@##.uhn-hpc.ca
```
where ## is your assigned student number. Again, make sure you replace BMB.pem with whatever the file name is that you chose!

It will then probably ask you if you're sure you want to connect. Type "yes" and press enter. Now you should see something like ```(base) ubuntu@ip-10-0-1-199:~$``` to the left of your cursor. Great, you're connected!

### Windows

First off, you'll need to save the .ppk key file that you've been given.

For Windows, we'll be logging in with Putty. So next, go to the "Session" category and fill in the "Host Name (or IP address)" field with:
```
ubuntu@##.uhn-hpc.ca
```
where ## is your assigned student number.

Now, in the "Connection" category, click the + next to "SSH". Click the + next to "Auth" and then click "Credentials". In the "Private key file for authentication" field, click browse and fine 

In the left hand categories, in the Connection category next to SSH click on the +. Click on + next to Auth and then click Credentials. In the private-key file for authentication field, hit browse and find the CBW.ppk file that you downloaded.

Now, click on the "Session" category again. In the "Saved sessions" field, type "Amazon node" and click Save. 

Now that Putty is set up, all you have to do is start putty and double-click on "Amazon node" to login. Now you should see something like ```(base) ubuntu@ip-10-0-1-199:~$``` to the left of your cursor. Great, you're connected!

## 1.2 Creating directories and moving around on the command line

Now that we're logged in, we're going to get familiar with what's on the server. First of all, type ```ls``` and press enter. You should see five things printed out:\
```CourseData  R  anaconda3  bin  workspace```\
These are all directories/folders on the server. Think of these like you would the directories on your own computer, e.g. "My Documents" and "Downloads". The ```ls``` command is just telling the server to list out the files/directories in the directory that we're currently in. 

Next, type in ```pwd``` and press enter. The ```pwd``` command tells you which directory you're currently in, so we can see that we're in ```/home/ubuntu```. 

Now we're going to change directory. We can change to any directory that exists here, but since we're mainly going to be working from ```workspace```, we'll change to there. Type in ```cd workspace``` and press enter (you will always need to press enter to run a command). If you type in ```pwd``` again, you should see that you've changed directory. If you type in ```ls```, you should see that the directory is empty. We're going to create directories for each module that we do, so we'll make one called ```bmb_module1```:\
```mkdir bmb_module1```\
Note that it doesn't really matter what we call this directory, but if you name it something different then you'll need to remember that when we run things later on so it's easier to keep it consistent.

Now let's change to the directory we just made:\
```cd bmb_module1```\
If you use the ```ls``` command, you should see that it's empty. Now what happens if we want to go back to the ```workspace``` directory? We can use ```cd ..```, which will take us up one level. Try that and then use the ```ls``` command again. You should see the folder that you made. If you use ```cd``` on it's own, then it'll take you back to where you started. Try it and then use the ```pwd``` command. Now to get back to the directory we created for this module, use: ```cd ~/workspace/bmb_module1```.
Note that any time you log out from the server (or lose connection/time out and get logged out), you will need to change back to the directory that you're working from. Use the ```pwd``` command once more to check that you're in the right place (```/home/ubuntu/workspace/bmb_module1```). 

## 1.3 Use wget to download files

Now we're going to download some files to use. There are a lot of different data repositories available online, but we're going to use some data from the Human Microbiome Project. This repository is quite straightforward because all we need is a link to the data and we can download it, but other repositories require their own programs and file formats for downloading which can make them quite complicated (and frustrating, if you don't know what you're doing!) to use. You can see the webpage that we're taking the files from [here](https://www.ibdmdb.org/downloads/).

Now we'll download the files. There are a few different commands that we can use for downloading files, so you might have seen others, or may see others in the future, but for now we'll be using a command called ```wget```. If you just type in ```wget``` on its own then you should see some information about it. You'll see that it gives an error message, because we haven't also given it a URL, but it also tells us that we can run ```wget --help``` to see more options. Try running that. If you scroll back up to where you ran it, you'll see a huge list of options that you could give ```wget``` depending on what you want to do. We won't be using most of these options for now, but you can usually add ```--help``` to whatever command you are trying to run to get some more information about what information it is expecting from you (these are called "arguments"). 

We're going to download three files like so:
```
wget https://www.ibdmdb.org/downloads/raw/HMP2/16S/2018-01-08/206534.fastq.gz
wget https://www.ibdmdb.org/downloads/raw/HMP2/16S/2018-01-08/206536.fastq.gz
wget https://www.ibdmdb.org/downloads/raw/HMP2/16S/2018-01-08/206538.fastq.gz
```
You should see some progress bars come up, but these files aren't very big so they shouldn't take very long to download. Now use ```ls``` again to see the files. You should see:
```
206534.fastq.gz  206536.fastq.gz  206538.fastq.gz
```

## 1.4 Move these files to the directories

When we downloaded these, we didn't make a directory to put them in, so let's do that now so that we can tidy them up a bit:
```
mkdir test_data
```
And then we can move them to this directory we've just made using the ```mv``` command. The ```mv``` command is expecting the name of the file that we want to move, and then the directory/path that we want to move this file to as arguments:
```
mv 206534.fastq.gz test_data/
```
If you use the ```ls``` command again now, you will see that there's only two files (along with the ```test_data``` directory). You can also use the ```ls``` command on the ```test_data``` directory: ```ls test_data```, and you should see the file that we moved into there.

Often, we might have a lot of files to move and we don't want to have to move them one by one. If we want to move multiple files of the same type, we can do that like this:
```
mv *.fastq.gz test_data/
```
The asterisk (```*```) acts as a wildcard, and anything that ends in ```.fastq.gz``` will be moved using this command. If you take a look in ```test_data``` again using the ```ls``` command, you should see all three files in there now, and not in your current directory. Another useful thing that you can do with the ```ls``` command is add another argument to get some information about the files: ```ls -lh test_data/```. When you run this, you should see who has permission to read/write to the files, the author and the owner of the files, the size of the files, and when they were created. The ```-l``` is what is telling this to give you the information, and adding the ```h``` converts this into "human" readable format. Try using it without the ```h``` - you'll see that the file sizes are in bytes instead of megabytes (MB/M).

## 1.5 Zip and unzip these files

You might have noticed the ```.gz``` on the end of the files indicating that they're zipped (or compressed). Sometimes when we run bioinformatics programs they can uncompress the files within the program, but other times we need to uncompress (unzip) them first. These are quite small files so you might think it's unnecessary to zip/compress them, but often we have thousands of sequencing data files and they can each be hundreds of gigabytes (GB) or even terabytes (TB) in size, so it becomes quite important to keep them compressed until we need them. 

We'll use the ```gunzip``` command to unzip them. Try typing in ```gunzip test_data/20``` and then pressing the tab button. If you press it a couple of times, you should see a list of your options come up. The tab button can be really useful for completing file paths for you and saving you a lot of typing! Continue typing and choose one of the files to unzip (e.g. type in ```8``` and then press tab again). 
Now if you run ```ls -lh test_data``` again, you should see that the file you unzipped no longer has the ```.gz``` on the end, and it's much larger in size now. 

We'll zip the file back up for now: ```gzip test_data/20``` - if you press tab to complete again, you should find that this will auto-fill with the file that you unzipped, because it's the only file type that the command is able to work with. 

What happens if you try to run this on a file that is already zipped? ```gzip test_data/206538.fastq.gz``` It should tell you that it's unable to run because ```.gz``` is already on that file. 

Let's unzip all of the files now: ```gunzip test_data/*```
See that we can use the asterisk (```*```) again as a wildcard and it will unzip every file in the ```test_data``` directory. Take a look at the directory with ```ls``` or ```ls -lh``` if you like.

## 1.6 Create new tar archive of files

There are several different ways of compressing files - there is ```gzip```/```gunzip``` that we just used, but we can also package up multiple files inside a directory together. We'll be using the ```tar``` command for this, and as you can see if you run ```tar --help```, you'll see that there are lots of options available for this. Let's try it out with the ```test_data``` directory:
```tar -czvf test_data.tar.gz test_data/```
Here we gave the command several arguments: ```-czvf``` (see below), ```test_data.tar.gz``` (the file name that we want our tar archive to have) and ```test_data/``` (the name of the directory that we want to compress. 
```-czvf``` is a short way of giving several arguments to the command: ```c``` for "create" (creates a new archive), ```z``` for "gZip" (this tells tar to write/read through gzip), ```v``` stands for "verbose" (meaning that it will print out information about what it is doing), and ```f``` for "file" or "archive". 

If you now run ```ls -lh```, you should see that the tar archive (```test_data.tar.gz```) is a smaller size than the 3 files would be together (check by running ```ls -lh test_data```). You can also take a look at what's in the tar archive with the ```less``` command: ```less test_data.tar.gz```. Press ```q``` (for "quit") when you're ready to exit. 

Usually we'll make a tar archive because we want to keep our files but save some space, so let's delete the original folder: ```rm -r test_data/```. Hopefully by now you're getting the hang of how these commands work. The ```rm``` command is for removing files - you should always be really careful when using it because it won't ask you if you're sure like your regular computer would, and most servers don't have a "recycle bin" for the files to be sent to, so if you remove them, they're gone for good. The ```-r``` argument is for removing a directory rather than just a file. 

## 1.7 Unzip tar archive

If we need to use the data that we zipped into the tar archive again, we'll need to unzip - or extract - it. 

To unzip the tar archive, we can do that like so: ```tar -xvf test_data.tar.gz```. Note that we just replaced the ```cz``` with ```x``` for "extract". You should be able to see the files back in ```test_data``` with the ```ls``` command.

## 1.8 Look at fasta and fastq files with less

Now we're going to take a look at these files. Let's look at ```206538``` first: ```less test_data/206538.fastq```. You can scroll through the file, and remember to press ```q``` when you want to stop looking at the file. If you want to look at it again, press the up arrow key to run the ```less``` command again. You can always press the up arrow to go back through the commands that you've run previously. If you've typed something wrong and want to start again, press ```ctrl```+```c``` to get back to a blank command prompt. 

You should have noticed that this file had the ```.fastq``` extension. Let's copy across the same files in ```fasta``` format: ```cp -r ~/CourseData/MIC_data/BMB_data/test_data_fasta/ .```
Here, we're using the ```cp``` command to copy across some files that I already set up earlier. The ```~``` is usually a shortcut for your home directory to save us typing it each time, and the rest is directories that we've already set up. The ```-r``` argument is the same as above for ```rm``` - it means that we're taking a directory and not a file, and then the ```.``` shows that we want to copy the data into the directory that we are currently in. We could replace it with another file path if we wanted. 

Take a look at the same file in ```fasta``` format: ```less test_data_fasta/206538.fasta```
You can also download these files by going to http://##.uhn-hpc.ca/ (where ## is your personal number). Navigate to the correct directory and then right click the files > save link as > choose location > choose name. You'll be able to open them with a text editor like TextEdit (Mac) or Notepad (Windows).
You should see that the fastq file has 4 lines for each sequence, while the fasta only has two. The fasta has a first line that starts with ">" that contains the sequence name and description, and the second line contains the actual DNA sequence (DNA in this case, but this file format could also contain RNA or protein sequences). fastq files start the same, with the first line containing the sequence name and description (but starting with an "@" symbol instead), and the second containing the sequence. The third line then contains a "+" character, and the fourth contains quality information about each base of the sequence (and should contain the same number of characters as the sequence). You can read more about the quality information that the symbols encode [here](https://en.wikipedia.org/wiki/FASTQ_format).

To count the number of lines in a file, we can use the ```less``` command again, with some additional information: ```less test_data_fasta/206538.fasta | wc -l``` and ```less test_data/206538.fastq | wc -l```

You should see that the fastq file contains double the number of lines as the fasta file. There are also other ways to count the number of sequences in a file, and these can be adapted for other purposes, too. E.g.: ```grep -c ">" test_data_fasta/206538.fasta``` - the ```grep``` command pulls out every occurrence of a phrase (or "string", as it's usually called in programming) and the ```-c``` argument tells it to count these. What happens if you don't use the ```-c``` argument? Why do you think this happened? 

In most programming languages, you have "positional" arguments and "named" arguments. Positional arguments need to be included in the proper position, or order. The order of positional arguments is defined within the program. Named or keyword arguments are given or passed to the program only after the name is given. In the case above, the ```-c ">"``` is a named argument and the file name ```test_data_fasta/206538.fasta``` is a positional argument.

**Question 1**: What happens if you try to do the same thing with "@" for the fastq file? Why is this? 

## 1.9 Installing programs to the server

Now we're going to learn how to install programs to the server. A lot of the commands we have just used (like ```grep``` and ```less```) are standard ones that will be installed on most servers, but frequently we will need to install other programs (usually called "packages"), too. The packages that we use are often not as stable as those that we download and use on our laptops (like Microsoft Word or Adobe Acrobat Reader) and so they sometimes depend on a particular version of another package. It frequently takes more time to install packages than it does to run them, and any bioinformatician will tell you how frustrating it can be. Anaconda can help to manage this, although it doesn't overcome these problems entirely! We already have it installed here, but you can follow the instructions [here](https://docs.anaconda.com/free/anaconda/install/linux/) to install it. This is how we installed it (NOTE: you don't need to run this part, but it won't hurt if you do!)
```
curl -O https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh
bash Anaconda3-2024.02-1-Linux-x86_64.sh
#hold down enter key until you get to the end of the agreement
#type yes
#confirm location by pressing enter
#yes
#now close and reopen the window - you'll need to log back in the same way as you did before!
```

## 1.10 Conda environments

Anaconda, or conda, allows us to have separate "environments" for installing packages into. This means that if one package requires version 3.1 of another package, but another requires version 2.9, they won't interfere with each other. Often when we're installing new packages or starting a new project, we'll make a new environment. This also helps us to keep track of which versions of a package we've used for a specific project. The environment is essentially a directory that contains a collection of packages that you've installed, so that other packages know where to access the package that they need to use. We're going to make a new environment to install some packages into:
```
conda create -n bmb_module1
```

You'll see that here we're using the ```conda``` command first, and then giving it the ```create``` and ```-n bmb_module1``` arguments. We could call this environment anything we like, but it's best to make this descriptive of what it is so that when we collaborate with others or share our code, it'll be obvious what this environment is for. 

Running this command will give you a warning:
```
WARNING: A conda environment already exists at '/home/ubuntu/CourseData/MIC_data/.conda/envs/bmb_module1'
Remove existing environment (y/[n])?
```
Type in ```y``` and press enter. 

You'll need to press ```y``` again at some point, to confirm that you want to install new packages. 

Now we can "activate" this environment like this: ```conda activate bmb_module1```.
Any time you are logged out and you log back in, you'll need to reactivate the environment if you want to be working from it. If you want to see the other environments that are installed and available, you can run ```conda info --envs```. You'll see that we've already installed a lot of environments that we'll be using over the next few days. While it would be great to be able to get you to install all of this for yourself, it can take quite a long time to install some packages so we've set most of them up already. If you want to see how we did that, you can see that [here](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2024:-Initial-setup-of-environments).

## 1.11 Install fastqc and multiqc

Now we'll install the packages that we want to use. Usually if there's a package that you're interested in, for example we'll be using one called "fastqc", you can just google "conda install fastqc" and you should see an anaconda.org page as one of the top hits, telling you how to install it. Sometimes you'll also see bioconda documentation, or a "package recipe" and this might give more details if you're struggling to install it. We'll install fastqc like this:
```
conda install bioconda::fastqc
```
You'll need to confirm that you want to install things with ```y``` at some point. If you forgot to activate the environment (see above), then you'll get an error that you don't have permissions to do this! 

You can test to see whether it got installed by typing ```which fastqc``` - this should show you the location that it's installed in. 

Now we'll install the second package that we need:
```
conda install bioconda::multiqc
```
Confirm this again with ```y```

As you might have guessed from the "qc" in both of these names, we'll be using them for Quality Control of the sequence data. 

## 1.12 Perform quality control on fastq files

First we'll be running fastqc, and to do that, we'll first make a directory for the output to go: ```mkdir fastqc_out```

Now we'll run fastqc:
```
fastqc -t 4 test_data/*.fastq -o fastqc_out
```
Here the arguments that we're giving ```fastqc``` are:
- ```-t 4```: the number of threads to use. Sometimes "threads" will be shown as ```--threads```, ```--cpus```, ```--processors```, ```--nproc```, or similar. Basically, developers of packages can call things whatever they like, but you can use the help documentation to see what options are available. We're using 4 here because that's the maximum that we have available. See below (```htop```) for how we find out about how many we have available.
- ```test_data/*.fastq```: the fastq files that we want to check the quality of.
- ```-o fastqc_out```: the folder to save the output to. 

## 1.13 htop - looking at the number of processes we have available or running

Try running ```htop```. This is an interactive viewer that shows you the processes that are running on your computer/server. There are a lot of different bits of information that this is showing us - you can see all of that [here](https://monovm.com/blog/what-is-htop-and-what-does-it-do/#:~:text=Htop%20is%20an%20interactive%20system,system%20processes%20can%20be%20viewed.), but the key things for us are:
- The CPUs (labelled 0, 1, 2, 3 at the top left) - this shows the percentage of the CPU being used for each core, and the number of cores shown here is the number of different processes/threads that we have available to us. In our case, this is 4. 
- Memory - this is the amount of memory, or RAM, that we have available to us. You'll see that it is ~16GB - this is similar to many laptops now, but many servers that you'll use or have access to for bioinformatics analysis will have much more than a standard computer. For example, one of our lab servers has ~1.5 TB RAM. The larger your dataset, or the deeper your sequencing depth, the more RAM you are likely to need.
- The processes (at the bottom) - you can see everything that is running under a PID (Process ID). This is useful when you're using a shared server to see who is running what, particularly for when you're wanting to run something that will use a lot of memory or will take a long time and you want to check that it won't bother anyone else. 

When you're done looking at this, press F10 (on a Mac this is fn+F10) to exit from this screen. 

## 1.14 Back to the quality control

Now take a look at one of the .html files in ```bmb_module1/fastqc_out/```

Next we'll run multiqc. The name suggests it might be performing QC on multiple files, but it's actually for combining the output together of multiple files, so we can run it like this:
```
multiqc fastqc_out --filename multiqc.html
```
So we've given as arguments:
- ```fastqc_out```: the folder that contains the fastqc output.
- ```--filename multiqc.html```: the file name to save the output as. 

Now look at multiqc.html.

There are some questions here to help you look at the files and interpret these:
- **Question 2**: What is the GC% of the samples?
- **Question 3**: What % of the samples are duplicate reads? Is this what you expected?
- **Question 4**: Now look at the Sequence Counts section. Which sample has the most reads?
- **Question 5**: How many unique and duplicate reads are in sample 206536?
- **Question 6**: Look at the Sequence Quality Histograms. Do these seem good to you? Why or why not? Does this seem normal?
- **Question 7**: Look at the top overrepresented sequence. If you want to see what it is, paste it into the "Enter accession number(s), gi(s), or FASTA sequence(s)" box [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and click on the blue "BLAST" button at the bottom of the page.\

## 1.15 Quality control on the data that we'll be using

Next, if you have time, you can try running the quality control on the data that we're going to be using for the other modules.

First of all, you can create a symlink to the data like this:
```
ln -s ~/CourseData/MIC_data/BMB_data/raw_data/ .
```
Creating a symlink is a bit like creating a desktop shortcut for something - we don't need to create a copy of it, but we just want easy access to it. You'll see that the format of the command is really similar to some of the copy (```cp```) or move (```mv```) commands that we used previously. If you use the ```ls``` command now, you should see that the directory that we've created the symlink to is a different colour to the other directories and files in the folder. These colours vary depending on the settings you choose in your console window, but they are usually different to indicate that this is a link. 

Now we'll create a directory for the fastqc output:
```
mkdir fastqc_raw_data_out
```

And run fastqc again:
```
fastqc -t 4 raw_data/*.fastq.gz -o fastqc_raw_data_out
```
Note that fastqc can run on zipped or unzipped files. 

And finally, run multiqc again:
```
multiqc fastqc_raw_data_out --filename multiqc_raw_data.html
```
Note that files with _R1 are forward reads and those with _R2 are reverse.
- **Question 8**: What do you notice between the R1 and R2 files? Why do you think that might be?
- **Question 9**: Do all of the samples have the same number of reads in the forward and reverse files?
- **Question 10**: When you look at the Sequence Quality Histograms you'll see two clusters. If you hover over the lines you can see the sample names. What do you notice?
- **Question 11**: Have a look again at the top overrepresented sequence. What is it?

## 1.16 Install QIIME2

For the final part of this module, we're going to install the program (QIIME2) that we'll be using for the rest of today and some of tomorrow. If you don't have enough time for this part, don't worry - we've also pre-installed this in a backup environment. 

First, remove the existing QIIME2 environment:
```
conda env remove -n qiime2-amplicon-2024.2
```

If you've found this easy so far and want more of a challenge, try doing this for yourself by following the instructions on the webpage [here](https://docs.qiime2.org/2024.2/install/).
**Hint**: You want to install the Amplicon version in a conda environment and the servers we're using are Linux.

Otherwise, follow the instructions here.

First, deactivate the conda environment that we're in now:
```
conda deactivate
```

Download the files that we'll need:
```
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml
```

Create an environment using this file:
```
conda env create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-linux-conda.yml
```
Note that where previously, we created an environment and then installed packages, this time we're creating an environment using a file that tells conda exactly how to create that environment and which packages should be installed within it. 

Once we've made it, we can activate this new environment:
```
conda activate qiime2-amplicon-2024.2
```

And now we can remove the file that we used to install it:
```
rm qiime2-amplicon-2024.2-py38-linux-conda.yml
```

## Answers

**Question 1**: What happens if you try to do the same thing with "@" for the fastq file? Why is this?\
The number of "@" in the fastq file is much more than the number of lines. This is because the "@" symbol is also used in the quality information. We can get round this by using part of the sample name, *e.g.,* ```grep -c "@206534" test_data/206534.fastq```.

**Question 2**: What is the GC% of the samples?\
51%

**Question 3**: What % of the samples are duplicate reads? Is this what you expected?\
In the "General Statistics" section, we can see that ~97% of the reads are duplicated. Looking in the "Sequence Counts" section and hovering over each sample will show us how many of the reads are unique. This makes sense, because the reads are from PCR-amplified samples so we are expecting most to occur more than once. 

**Question 4**: Now look at the Sequence Counts section. Which sample has the most reads?\
206354.

**Question 5**: How many unique and duplicate reads are in sample 206536?\
752 and 29,883.

**Question 6**: Look at the Sequence Quality Histograms. Do these seem good to you? Why or why not? Does this seem normal?\
The sequence quality here is all really high! These are all good sequences, but this isn't normal. This is because the samples that are available for download from the HMP website have already been quality filtered. 

**Question 7**: Look at the top overrepresented sequence. If you want to see what it is, paste it into the "Enter accession number(s), gi(s), or FASTA sequence(s)" box [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and click on the blue "BLAST" button at the bottom of the page. \
All of the top hits are Bacteroides (finegoldii, caccae, stercoris, etc). These samples are from HMP2 IBD gut samples, so this seems normal!

**Question 8**: What do you notice between the R1 and R2 files? Why do you think that might be?\
There are higher numbers of duplicated reads in the reverse (R2) read files. 

**Question 9**: Do all of the samples have the same number of reads in the forward and reverse files?\
Yes! This is good, because they should all be paired reads. If they were different then this would indicate that something in our sequencing had gone wrong.

**Question 10**: When you look at the Sequence Quality Histograms you'll see two clusters. If you hover over the lines you can see the sample names. What do you notice?\
The group of lines with lower mean quality scores corresponds to the reverse (R2) reads. The quality for reverse reads from a MiSeq sequencing run is always a little lower than for the forward reads.

**Question 11**: Have a look again at the top overrepresented sequence. What is it?\
*Cystobacter* sp. is the top hit, and *Alsobacteraceae* bacterium R-9 is the second hit.
