#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/Geoduck_TagSeq/output/Stringtie2

#load packages
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.1-GCCcore-7.3.0 #Transcript assembly: StringTie
module load gffcompare/0.11.5-foss-2018b #Transcript assembly QC: GFFCompare

stringtie --merge -p 8 -G ../../../refs/Panopea-generosa-v1.0.a4.mRNA_SJG.gff3 -o Pgenerosa_merged.gtf gtf_list.txt #Merge GTFs to form $
echo "Stringtie merge complete" $(date)

gffcompare -r ../../../refs/Panopea-generosa-v1.0.a4.mRNA_SJG.gff3 -G -o merged Pgenerosa_merged.gtf #Compute the accuracy and pre$
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

python ../../scripts/prepDE.py -g Pgenerosa_gene_count_matrix.csv -i listGTF.txt #Compile the gene count matrix
echo "Gene count matrix compiled." $(date)

