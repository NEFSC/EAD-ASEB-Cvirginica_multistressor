#!/bin/bash
#SBATCH --job-name="hisat2_align"
#SBATCH -t 048:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"

# before running..
# make directory named hisat2 as output/hisat2/

cd #nav back to home directory (allows job to be run from anywhere)

# load modules, requires hisat2 and samtools
module load bio/hisat2/2.2.1
module load bio/samtools/1.11

cd Cvirginica_multistressor_TagSeq/output/hisat2 # nav to hisat2 - symbolic directory works well when output to the current dir as ./

# symbolically link clean reads to hisat2 dir
ln -s ../fastp_multiQC/clean/*.fastq.gz ./ # call the .fastq.gz output from fastp trim - make symb link to output/hisat2

# activate python fro hisat2-build
source ../../../python_venv/bin/activate  # activate python cirtual envriomment to call python and run hisat2-build

# index the reference genome for Panopea generosa output index to working directory
hisat2-build -f ../../../refs/GCF_002022765.2_C_virginica-3.0_genomic.fna ./Cvirginica_ref
echo "Referece genome indexed. Starting alingment" $(date)

# exit python virtual envrionment
deactivate

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        hisat2 -p 8 --dta -x Cvirginica_ref -U ${i} -S ${i}.sam
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
                echo "${i} bam-ified!"
        rm ${i}.sam
done
