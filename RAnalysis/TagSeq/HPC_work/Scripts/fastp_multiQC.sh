#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/Geoduck_TagSeq/output/fastp_multiQC

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-2.7.15

# symbolically link 'clean' reads to hisat2 dir
ln -s ../../../../KITT/hputnam/20201217_Geoduck_TagSeq/*.fastq.gz ./

# Make an array of sequences to trim
array1=($(ls *.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 
        fastqc clean.${i}
done 

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

multiqc ./ #Compile MultiQC report from FastQC files

echo "Cleaned MultiQC report generated." $(date)
