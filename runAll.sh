#!/bin/bash
#######################
# VirusSeqPipeline
# By Roseric Azondekon
#######################

# set working directory
folder=$(pwd) #capture the VirusSeq working directory
cd "$folder"


# Get Reference genome files if they do not exist
if [ ! -f "$folder"/Mosaik_JumpDb/hg19.fa ] ; then
echo "Downloading hg19.fa.gz..."
wget "$folder" "http://odin.mdacc.tmc.edu/~xsu1/hg19.fa.gz"
echo "Extracting hg19.fa.gz"
gunzip -c hg19.fa.gz > "$folder"/Mosaik_JumpDb/hg19.fa
rm -f hg19.fa.gz
fi

if [ ! -f "$folder"/Mosaik_JumpDb/gibVirus.fa ] ; then
echo "Downloading gibVirus.fa.gz..."
wget "$folder" "http://odin.mdacc.tmc.edu/~xsu1/gibVirus.fa.gz"
echo "Extracting gibVirus.fa.gz"
gunzip -c gibVirus.fa.gz > "$folder"/Mosaik_JumpDb/gibVirus.fa
rm -f gibVirus.fa.gz
fi

if [ ! -f "$folder"/Mosaik_JumpDb/hg19Virus.fa ] ; then
echo "Downloading hg19Virus.fa.gz..."
wget "$folder" "http://odin.mdacc.tmc.edu/~xsu1/hg19Virus.fa.gz"
echo "Extracting hg19Virus.fa.gz"
gunzip -c hg19Virus.fa.gz > "$folder"/Mosaik_JumpDb/hg19Virus.fa
rm -f hg19Virus.fa.gz
fi

#Check if all Jump files exist. If not, create them
exts=( "$folder"/Mosaik_JumpDb/*.jmp ); 
JMP=${#exts[@]}
if [ "$JMP" -eq 9 ]; then
echo "All JMP files already exist...Skip Jump Db creation..."
else
cd Mosaik_JumpDb
bash Jump_file_builder.sh
fi


#For each sample in the folder samples
folder2="$folder"/samples
cd "$folder2"
for sample in *; do
# All job processing goes here
##create folder if not exist
mkdir -p "$folder2"/"$sample"/Gig
mkdir -p "$folder2"/"$sample"/SV_gDNA
mkdir -p "$folder"/results/"$sample"



##################################################################################
#1. Virus detection by NGS data
##################################################################################

##Converting external read formats to the native MOSAIK format
"$folder"/Mosaik_bin/MosaikBuild -q "$folder2"/"$sample"/"$sample"_1.fq.gz -q2 "$folder2"/"$sample"/"$sample"_2.fq.gz -out "$folder"/results/"$sample"/"$sample"_Virus.bin -st illumina

##performing alignment against human genome reference hg19 with MosaikAligner
"$folder"/Mosaik_bin/MosaikAligner -in "$folder"/results/"$sample"/"$sample"_Virus.bin -ia "$folder"/Mosaik_JumpDb/hg19.fa.bin -out "$folder"/results/"$sample"/"$sample"_Virus.bin.aligned -hs 15 -mmp 0.1 -mmal 0.5 -act 25 -mhp 100 unique -j "$folder"/Mosaik_JumpDb/hg19.JumpDb -p 14 -km -pm -rur "$folder"/results/"$sample"/"$sample"_Unalg.fq

#start aligning unmapped reads against virus genomes to detect the virus

##sorting the mapped reads in terms of their mapped genomic location using MosaikSort
"$folder"/Mosaik_bin/MosaikBuild -q "$folder"/results/"$sample"/"$sample"_Unalg.fq -out "$folder"/results/"$sample"/"$sample"_Virus.bin -st illumina

"$folder"/Mosaik_bin/MosaikAligner -in "$folder"/results/"$sample"/"$sample"_Virus.bin -ia "$folder"/Mosaik_JumpDb/gibVirus.fa.bin -out "$folder"/results/"$sample"/"$sample"_Virus.bin.aligned -hs 15 -mmp 0.15 -act 25 -mhp 100 -m all -j "$folder"/Mosaik_JumpDb/gibVirus.JumpDb -p 14 -km -pm

"$folder"/Mosaik_bin/MosaikSort -in "$folder"/results/"$sample"/"$sample"_Virus.bin.aligned -out "$folder"/results/"$sample"/"$sample"_Virus.bin.aligned.sorted

##producing an assembly of the reads pileup for visualization purposes using MosaikAssembler
"$folder"/Mosaik_bin/MosaikAssembler -in "$folder"/results/"$sample"/"$sample"_Virus.bin.aligned.sorted -ia "$folder"/Mosaik_JumpDb/gibVirus.fa.bin -out "$folder2"/"$sample"/Gig/"$sample"_Virus.bin.aligned.sorted.assembled -f ace > "$folder"/results/"$sample"/"$sample"_VirusLog.txt

#detect the virus in the sample by VirusSeq_Detection.pl
perl "$folder"/VirusSeq_Script/VirusSeq_Detection.pl "$folder"/results/"$sample"/"$sample"_VirusLog.txt 1000 "$folder"/results/"$sample"/"$sample"_VirusName.txt



##################################################################################
#2. Detection of virus integration sites by NGS data
##################################################################################

##################################################################################
#2.1 PE reads FASTQ file mapping against hybrid genome hg19Virus
##################################################################################

#start alignment against hg19Virus
"$folder"/Mosaik_bin/MosaikBuild -q "$folder2"/"$sample"/"$sample"_1.fq.gz -q2 "$folder2"/"$sample"/"$sample"_2.fq.gz -out "$folder"/results/"$sample"/"$sample".bin -st illumina

"$folder"/Mosaik_bin/MosaikAligner -in "$folder"/results/"$sample"/"$sample".bin -ia "$folder"/Mosaik_JumpDb/hg19Virus.fa.bin -out "$folder"/results/"$sample"/"$sample".bin.aligned -hs 15 -mmp 0.06 -mmal -minp 0.5 -act 25 ‐mhp 100 -m unique -a all -j "$folder"/Mosaik_JumpDb/hg19Virus.JumpDb -km -pm -p 14


#Spanner is used to generate the list of discordant reads across different chromosomes
"$folder"/Mosaik_bin/Spanner --scan --infile "$folder"/results/"$sample"/"$sample".bin.aligned --outdir "$folder2"/"$sample"/SV_gDNA

cd "$folder2"/"$sample"/SV_gDNA

for i in ./*"$sample".bin*;do mv -- "$i" "${i//"$sample".bin/MSK}";done

cd "$folder"

"$folder"/Mosaik_bin/Spanner --build --infile "$folder"/results/"$sample"/"$sample".bin.aligned --outdir "$folder2"/"$sample"/SV_gDNA -f "$folder2"/"$sample"/SV_gDNA/MSK.stats -a "$folder"/Mosaik_JumpDb/Spanner_anchor_hg19Virus.txt -t

##################################################################################
#3.2 Single‐end reads FASTQ file mapping against hybrid genome hg19Virus
##################################################################################

#align the first mate against hg19Virus
"$folder"/Mosaik_bin/MosaikBuild -q "$folder2"/"$sample"/"$sample"_1.fq.gz -out "$folder"/results/"$sample"/"$sample"_1.bin -st illumina

"$folder"/Mosaik_bin/MosaikAligner -in "$folder"/results/"$sample"/"$sample"_1.bin -ia "$folder"/Mosaik_JumpDb/hg19Virus.fa.bin -out "$folder"/results/"$sample"/"$sample"_1.bin.aligned -hs 15 -mmp 0.06 -mmal -minp 0.5 -act 25 -mhp 100 -m unique -j "$folder"/Mosaik_JumpDb/hg19Virus.JumpDb -p 14 -km -pm

"$folder"/Mosaik_bin/MosaikSort -in "$folder"/results/"$sample"/"$sample"_1.bin.aligned -out "$folder"/results/"$sample"/"$sample"_1.bin.aligned.sorted -u

#align the second mate against hg19Virus
"$folder"/Mosaik_bin/MosaikBuild -q "$folder2"/"$sample"/"$sample"_2.fq.gz -out "$folder"/results/"$sample"/"$sample"_2.bin -st illumina

"$folder"/Mosaik_bin/MosaikAligner -in "$folder"/results/"$sample"/"$sample"_2.bin -ia "$folder"/Mosaik_JumpDb/hg19Virus.fa.bin -out "$folder"/results/"$sample"/"$sample"_2.bin.aligned -hs 15 -mmp 0.06 -mmal -minp 0.5 -act 25 -mhp 100 -m unique -j "$folder"/Mosaik_JumpDb/hg19Virus.JumpDb -p 14 -km -pm

"$folder"/Mosaik_bin/MosaikSort -in "$folder"/results/"$sample"/"$sample"_2.bin.aligned -out "$folder"/results/"$sample"/"$sample"_2.bin.aligned.sorted -u

#merge two sorted files
"$folder"/Mosaik_bin/MosaikMerge -in "$folder"/results/"$sample"/"$sample"_1.bin.aligned.sorted -in "$folder"/results/"$sample"/"$sample"_2.bin.aligned.sorted -out "$folder"/results/"$sample"/"$sample"_SE.bin.aligned.sorted

"$folder"/Mosaik_bin/MosaikText -in "$folder"/results/"$sample"/"$sample"_SE.bin.aligned.sorted -axt "$folder"/results/"$sample"/"$sample"_SE.bin.aligned.sorted.axt



##################################################################################
#3.3 Detection of virus integration site
##################################################################################

#Spanner cross reads converter
cd "$folder2"/"$sample"/SV_gDNA

for i in ./*"$sample".bin*;do mv -- "$i" "${i//"$sample".bin/MSK}";done

perl "$folder"/VirusSeq_Script/Spanner_cross_converter.pl "$folder"/Mosaik_JumpDb/hg19Virus_refGene_RIS.txt "$folder"/results/"$sample"/"$sample"_SE.bin.aligned.sorted.axt "$folder"/results/"$sample"/"$sample"_CrossReads.txt 

##Virus integration site detection 
perl "$folder"/VirusSeq_Script/VirusSeq_Integration.pl "$folder"/results/"$sample"/"$sample"_CrossReads.txt "$folder"/Mosaik_JumpDb/hg19Virus_refGene_RIS.txt hg19 192 95 50 "$folder"/results/"$sample"/"$sample"_Integration.txt
done

################################### END ##########################################
