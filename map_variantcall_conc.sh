#! /bin/bash

######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to map RNAseq data to a set of CDS, call variants an take concensus 
##          sequence with reference replaced with new variants
##      Start with cleaned reads
## For running the script on the Alabama Super Computer.
##	For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## 	then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
##  You may need to increase these for bigger datasets
## 		queue: medium
##		core: 6
##		time limit (HH:MM:SS): 18:00:00 
##		Memory: 12gb
##		
###############################################


########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1
module load multiqc
module load trimmomatic/0.39
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread
#module load gffcompare

#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x


##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubtss          ## Example: MyID=aubtss


WD=/scratch/$MyID/PracticeRNAseq_Full_Script            ## Example:/scratch/$MyID/PracticeRNAseq  
CD=$WD/CleanData            				## Example:/scratch/$MyID/PracticeRNAseq/CleanData   #   *** This is where the cleaned paired files are located from the last script
REFD=$WD/Ref          ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome    # this directory contains the indexed reference genome for the garter snake
MAPD=$WD/Map_HiSat2           			## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2      #
RESULTSD=/home/$MyID/PracticeRNAseq_Full/Counts_H_S_2024      ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S
REF=IIS_CDS                 ## This is what the "easy name" will be for the genome reference
VAR=$WD/Variants

############################***********  Mapping and Calling SNPs ************##########################

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $VAR
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD
### Copy the reference set of cds (.fasta) to this REFD directory
cp /home/${MyID}/class_shared/Dog_Tasha_GCF_000002285.5/data/GCF_000002285.5/${REF}.fasta .


#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build ${REF}.fasta IIS_CDS_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
cd ${CD}  #### This is where our clean paired reads are located.

## Create list of fastq files to map.    Example file format of your cleaned reads file names: SRR629651_1_paired.fastq SRR629651_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: SRR629651

## Move to the directory for mapping
cd ${MAPD}

## move the list of unique ids from the original files to map
mv ${CD}/list  . 

## process the samples in the list, one by one using a while loop
while read i;
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x "${REFD}"/IIS_CDS_index       \
    -1 "${CD}"/"$i"_1_paired.fastq  -2 "${CD}"/"$i"_2_paired.fastq      \
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam  

    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 6  "$i".bam  -o  "$i"_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  	## remove unmapped reads
  samtools view -F 0x04 -b "$i"_sorted.bam > "$i"_sorted_mapOnly.bam


##############  Calling SNPs  #######################
   #Call SNPs
bcftools mpileup -f reference_genome.fa "$i"_sorted.bam | bcftools call -mv -Ov -o "$i"_variants.vcf

 #Generate the consensus sequence
bcftools consensus -f reference_genome.fa -o "$i"_IISconsensus.fasta "$i"_variants.vcf

done<list


#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *_variants.vcf ${RESULTSD}

### copy the final results files (the count matricies that are .cvs) to your home directory. 
cp *_consensus.fasta ${RESULTSD}
