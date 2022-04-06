#!/bin/bash
#SBATCH --job-name="crgKEGG_diamond"
#SBATCH -t 004:00:00
#SBATCH --mem=32GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output=./"%x_out.%j"
#SBATCH --error=./"%x_err.%j"

# job run (in SEDNA) from within the folder 'output/diamond_KEGGIDs' as the following
# sbatch ../../scripts/crgKEGG_diamond.sh ../../../refs/Cgigas_KEGG_IDs_translated.fasta ../../../refs/GCF_002022765.2_C_virginica-3.0_rna.fna

# objective: of this code is to obtain a list of blast hits of the Cvirginica genome to the C gigas protein database 
# this is specifically to obtain KEGG Ids since only Cgigas (currently - April 2022) has this annotation level 
# the prodct of this code will provide rationale (best blast hit, i.e. E-score) to merge KEGG Ids of Cgagis to C virginica
 

 
date # call curernt data, writes to output file to gauge the speed of diamond (proposed to be 100x faster than blast)

# Note: diamond is an alternative and more efficient module for blastx
module load bio/blast/2.11.0+
module load bio/diamond/2.0.14.152

# Pass the database you want as the first argument
database=$1 # the first argument called after the slurm .sh call (Ex: bash.sh filename1.fasta filename2.fasta - inthis case $1 is filename1.fasta)
# 1$ is the protein fasta file of KEGGIDs for Cgigas -  ../../../refs/KEGG_Crassostrea_gigas.fasta ('output/diamond_KEGGIDs' directory where sbatch from)

# Pass the query you want as the second argument
query=$2 # the second argumnet, review what this mean above 
# 2$ is the C virginica annotated mRNA file  - ../../../refs/GCF_002022765.2_C_virginica-3.0_rna.fna ('output/diamond_KEGGIDs' directory where sbatch from)

mkdir ./Cgigasdb/

# Make the database - note build a protein database and follow with blastx (nucl query w/ protein database)
diamond makedb --in $1 -d ./Cgigasdb/Cgigas_db

# runs blast on the P generosa gene fasta against the Pacific oyster protein database we created above
# --very sensitive finds hits with best sensitivity <40% identity 

diamond blastx -d ./Cgigasdb/Cgigas_db.dmnd -q $2 -o ./crgKEGG_diamond_out --outfmt 6 

# qseqid sseqid pident evalue length qlen slen qstart qend sstart send sseq

echo "Done"
date

