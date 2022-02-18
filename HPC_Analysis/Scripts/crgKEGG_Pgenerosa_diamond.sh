#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/refs

date

# Note: diamond is an alternative and more efficient module for blastx
# module load BLAST+/2.7.1-foss-2018a
module load DIAMOND/2.0.0-GCC-8.3.0

# Pass the database you want as the first argument
# Pass the query you want as the second argument

database=$1
query=$2

mkdir ./Cgigasdb/

# Make the database - note build a protein database and follow with blastx (nucl query w/ protein database)
diamond makedb --in $1 -d ./Cgigasdb/Cgigas_db


# runs blast on the P generosa gene fasta against the Pacific oyster protein database we created above
# --very sensitive finds hits with best sensitivity <40% identity 

diamond blastx -d ./Cgigasdb/Cgigas_db.dmnd -q $2 -o ./crgKEGG_diamond_out --outfmt 6 

# qseqid sseqid pident evalue length qlen slen qstart qend sstart send sseq

echo "Done"
date

