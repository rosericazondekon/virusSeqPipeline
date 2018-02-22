VirusSeqPipeline: A VirusSeq pipeline for processing multiple RNAseq samples: a command line approach 
=========

The VirusSeq software is a next generation sequencing software of human cancer tissue. Its goal is to identify viruses and their integration sites in human cancer cells. VirusSeq has originally been developed by [Chen and collaborators](https://academic.oup.com/bioinformatics/article/29/2/266/202055 "VirusSeq: software to identify viruses and their integration sites using next-generation sequencing of human cancer tissue") and reported in a 2012 paper published in Bioinformatics. Although the software is open source and freely available, its installation is not at the ease of the novice user. VirusSeq has not been designed to effectively and efficiently process multiple samples. In its present shape, it requires the user to write a shell script to every single sample to process. Furthermore, the user has to engage in the tedious and error prone tasks of path modification. All these limitations of the VirusSeq software may tremendously limit its use in the scientific community. 

In this project, we aim at addressing the aforementioned limitations by designing an easy to install, easy to implement pipeline for the VirusSeq software. Our approach does not require the user to have knowledge of any shell command. In this user guide documentation, we present the file organization and a more compact package for the user to install on his computer. 



VirusSeqPipeline description 
=========
This pipeline is composed of the MOSAIK aligner, reference genome and annotation files and PERL programs. The PERL programs are originally from the VirusSeq software. In the full version of this pipeline, we provided the user with the reference genome and annotation files made available by the authors of VirusSeq and available at : http://odin.mdacc.tmc.edu/~xsu1/VirusSeq.html. Also, see the original user guide at the same URL for more information. 

The file organization of our pipeline is presented as follow: 

- the directory VirusSeq/Mosaik_bin contains all the Mosaik aligner executable tools 

- the directory VirusSeq/MosaikJumpDb contains the following files: 

- the reference genome and annotation files: 

    1. hg19.fa  

    2. gibVirus.fa  

    3. hg19Virus.fa  

    4. hg19Virus_refGene_RIS.txt  

    5. Spanner_anchor_hg19Virus.txt 

- a SHELL script Jump_file_builder.sh to generate the jump files from the reference genome and annotation files 

- the VirusSeq/VirusSeq_Script directory contains the following PERL programs: 

    1. Spanner_cross_converter.pl 

    2. VirusSeq_Detection.pl 

    3. VirusSeq_Integration.pl 

- the VirusSeq/samples directory is where the user is supposed to place the samples to be processed. In section 2, we provide detailed direction on how the user is expected to organize his FASTQ samples. 

- the VirusSeq/results directory contains the files resulting from the processing of each sample organized in folders named after each sample. 

- The VirusSeq folder contains the file runAll.sh which is a SHELL script containing the command to run the VirusSeq pipeline on all the samples contained in the VirusSeq/samples directory. This directory also contain other archive files for the MOSAIK aligner, the reference genome and annotations file as well as a copy of the original VirusSeq package. In addition, we also provide a copy of the eagleView software used to preview the viruses sequences in .ACE format. On some computers and for users with limited privileges, the Mosaik tools might not be allow to run as executable programs. This usually happens when users attempt to run this pipeline from an external media. Knowing that this may limit the use of this pipeline, we provide the file AllowExec.txt to help users figure out how to allow the Mosaik toole to run as executable. 


We bring to the user's attention to the fact that altering or changing this file organization may prevent VirusSeq from running correctly. 



VirusSeqPipeline data input
=========
This pipeline requires the user to input RNA-seq data in a specific compressed FASTQ format (.fq.gz) with paired-end reads. If for example the user is to process a sample "L526401A", the FASTQ files "L526401A_1.fq.qz" and "L526401A_2.fq.gz" are required. We expect the user to create a folder named "L526401A" containing the files "L526401A_1.fq.qz" and "L526401A_2.fq.gz". The folder "L526401A" is later to be placed in the directory VirusSeq/samples for its processing. 

In general, Each folder in the VirusSeq/samples directory is a sample and a sample is made of two FASTQ files with the respective extensions "_1.fq.gz" and "_2.fq.gz". The folder in the VirusSeq/samples directory must be named after the name of the sample it represents. 

It is important to note that after the processing of each sample, two additional folders named respectively "Gig" and "SV_gDNA" are created inside in subfolder of the VirusSeq/samples directory. The "Gig" folder contains the sequence files (in .ACE format) of all the viruses reference found for that sample. Each .ACE file can be preview using the freely availabel EagleView package we provided along with this pipeline. EagleView is also available for download at: http://www.niehs.nih.gov/research/resources/software/biostatistics/eagleview/index.cfm. 

The "SV_gDNA" folder received the results from the virus integration sites. 



How to install and run VirusSeqPipeline
=========
Please, follow the instructions below to set up the pipeline and run it: 

- Clone VirusSeqPipeline

	git clone https://github.com/rosericazondekon/VirusSeqPipeline.git

- Now open a shell terminal (CTRL+TAB+T on Ubuntu) and change your directory (using the "CD" command) to your VirusSeqPipeline directory 

- Make sure that all the Mosaik aligner tools are allowed to be executed as files. To do this, navigate to the VirusSeq/Mosaik_bin directory and right click on each file, then select "Properties", next select "Permissions" and tick the box "Allow executing file as program". If for some reasons, this is impossible, follow the instructions in "AllowExec.txt" file. 

- Put your samples organized as explained in section 2 in the VirusSeq/samples directory. Do remember that each sample directory must contain two compressed FASTQ files ending respectively by the extensions "_1.fq.gz" and "_2.fq.gz" 

- And run the following shell command:
```shell
bash runAll.sh
```

Do notice that depending on the number of samples, you might want to make sure that hundreds of Gigabytes of space are available on your media. The entire processing may therefore take several hours or days to complete.




~
=========

Roseric Azondekon