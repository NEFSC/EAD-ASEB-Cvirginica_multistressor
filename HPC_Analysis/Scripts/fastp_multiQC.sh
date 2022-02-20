#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /sgurr/Cvirginica_multistressor_TagSeq/output/fastp_multiQC

# load modules needed
module load bio/fastp/0.23.2
module load bio/fastqc/0.11.9 

# symbolically link 'clean' reads to hisat2 dir
ln -s ../../../../../share/nefsc/mcfarland_sequecenes/TagSeq_oysters_2021/SA21200*/*.fastq.gz ./ # call backward from the directory to the share folder, input symbolic link to directory (-D in SBATCH header)

# Make an array of sequences to trim
array1=($(ls SA21200*/*.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 
        fastqc clean.${i}
done 

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

source ../../../python_venv/bin/activate # from the current directory (-D in SBATCH header), activates the bin of installed python packages, including multiqc

multiqc ./ #Compile MultiQC report from FastQC files

echo "Cleaned MultiQC report generated." $(date)
