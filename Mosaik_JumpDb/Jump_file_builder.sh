################################################
###Please use "qsub jma.sh" to submit to cluster
###############################################

#!/bin/bash
###
#PBS -S /bin/bash
#PBS -N JumpDb_builder
#PBS -d /RIS/home/xsu1/VirusSeq/Mosaik_JumpDb
#PBS -e /RIS/home/xsu1/VirusSeq/Mosaik_JumpDb/errfile.txt
#PBS -o /RIS/home/xsu1/VirusSeq/Mosaik_JumpDb/outfile.txt
###PBS -M xsu1@mdanderson.org 
#PBS -q batch
#PBS -l nodes=1:ppn=1

#####################
##your command here##
#####################

#set working directory
folder=$(pwd)
dir=${folder%/Mosaik_JumpDb}
dir="$dir"/Mosaik_bin

##for gib virus reference genome
$dir/MosaikBuild -fr gibVirus.fa -oa gibVirus.fa.bin -st illumina -assignQual 40
$dir/MosaikJump -ia gibVirus.fa.bin -out gibVirus.JumpDb -hs 15 -mhp 100 -km

##version-0.89 and for combined hg19Virus reference genome
$dir/MosaikBuild -fr hg19Virus.fa -oa hg19Virus.fa.bin -st illumina -assignQual 40
$dir/MosaikJump -ia hg19Virus.fa.bin -out hg19Virus.JumpDb -hs 15 -mhp 100 -km

##for hg19 reference genome
$dir/MosaikBuild -fr genome_ref.fa -oa genome_ref.fa.bin -st illumina -assignQual 40
$dir/MosaikJump -ia genome_ref.fa.bin -out genome_ref.JumpDb -hs 15 -mhp 100 -km
