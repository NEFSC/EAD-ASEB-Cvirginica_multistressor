# Eastern oyster (*C. virginica*) TagSeq Bioinformatics

## <span style="color:blue">**Table of Contents**</span>
  - [Upon upload to HPC...](#Initial-diagnostics-upon-sequence-upload-to-HPC)
	  - [Upon upload to HPC...](#Initial-diagnostics-upon-sequence-upload-to-HPC)
      - [Count raw reads](#Count-the-number-of-read-files)
      - [Digital fingerprint md5sums](#run-checksum) (HPC script; <span style="color:green">**md5_checksum.sh**<span>)
  - [1. MultiQC: Initial QC](#Quality-check-of-raw-reads) (HPC script; <span style="color:green">**mutliqc.sh**<span>)
  - [2. Trimming and QC of 'clean' reads](#Trimming-and-post-trim-quality-check-of-'clean'-reads)
  	- [Remember! polyA tail in TagSeq](#Trimming-polyA-tail)
	- [fastp - about/commands](#What-this-script-will-do...)
	- [fastp and MulitiQC: Trim and QC](#shell-script-fastp_multiqc.sh) (HPC script; <span style="color:green">**fastp_mutliqc.sh**<span>)
  - [3. Alignment of cleaned reads to reference](#HISAT2-Alignment-of-cleaned-reads-to-reference)
	- [Upload reference genome](#Reference-genome-upload-to-HPC)
	- [HISAT2 - about/commands](#HISAT2-alignment)
	- [samtools - about/commands](#samtools)
	- [HISAT2: Index reference and alignment](#HPC-Job-HISAT2-Index-Reference-and-Alignment) (HPC script; <span style="color:green">**HISAT2.sh**<span>)
  - [4. Assembly and quantification](#Assembly-and-quantification)
  	- [Upload annotation reference for assembly](#Upload-annotation-reference-gff-or-gff3-to-HPC)
	- [StringTie2 - about/commands](#StringTie)
	- [gffcomapare - about/commands](#gffcompare)
	- [prepDE.py - about/commands](#Python-step-prepDE.py) (essential prep, load open-source python script for count matrix)
		- [List HISAT2 output files (for --merge in Stringtie2)](#gtf_list.txt-run) (essential Stringtie2 ```--merge``` prep; **gtf_list.txt** file!)
		- [List HISAT2 output files (for prepDE.py)](#listGTF.txt-run) (essential prepDE.py prep; **listGTF.txt** file!)
	- [Stringtie2: Assembly step](#HPC-job-Assembly) (HPC script; <span style="color:green">**Stringtie2.sh**<span>)
	- [Stringtie2, gffcomapre, prepDE.py: Merge and build read count matrix for DEG analysis](#HPC-job-Merge-and-Build-Read-Count-Matrix-for-DEG-analysis) (HPC script; <span style="color:green">**Stringtie2_merge_prepDEpy.sh**<span>)


# Initial diagnostics upon sequence upload to HPC
--------------------------------------------
### Count the number of read files
ls -1 | wc -l

should equal 141 seq samples *2lanes per sample + 1 md5 = 283

**NOTE:** H.Putnam ran transfer_checks.sh and output into the 20201217_Geoduck_TagSeq/ folder


# run checksum

```
nano transfer_checks.sh
```
```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

# generate md5
md5sum /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/*.gz > /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/URI.md5

# Count the number of sequences per file

zcat /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/*fastq.gz | echo $((`wc -l`/4)) > /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/rawread.counts.txt
```

```
sbatch transfer_checks.sh
```

- check the digital fingerprint of the files with md5sum
- compare md5sum of our output URI.md5 file to the UT Library Information pdf; okay the upload - move forward






# Quality check of raw reads
-------------------------------------------





## fastqc

**About:** A quality control tool for high throughput sequence data. [website here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

outputs html report for quality check


## MultiQC

**About:** MutliQC aggregates the output html from fastqc into a cumulative interface to easily quality check/diagnose your sequence data pre/post trim! [website here](https://multiqc.info/)

Note! **multiQC is installed to sedna using python**, see calls  below for how to start a python virtual environment, install packages, and navigate to them  within slurm jobs

## Python virtual environment

start your personnel python python_venv
```
python3 -m venv /<your desired path to venv>/
```
after complete, activate the virtual environment with the following..
```
source <venv>/bin/activate
```
you can now download packages/modules that are otherwise unavailable such as multiQC.
use 'pip install' to download
```
pip install multiQC
```
you now have mutliqc installed in the <venv>/ bin folder
**to call python modules in slurm jobs** simply enter the virtual environment using 'source' <venv>/bin/activate and all modules can then be called. *you do not need to manually load modules installed to the python venv*

to leave the virtual envrionment...
```
deactivate
```


# shell script: <span style="color:green">**raw_multiQC,sh**<span>
```
#!/bin/bash
#SBATCH --job-name="multiQC_raw_reads"
#SBATCH -t 002:00:00
#SBATCH --mem=32GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"

# load modules needed
module load bio/fastp/0.23.2
module load bio/fastqc/0.11.9

# run sbatch from the sgurr home directory..

cd Cvirginica_multistressor_TagSeq/output/fastp_multiQC # run from the sgurr directory, navigate to the fast_multiQDC folder
mkdir raw # make directory for trimmed fastq files and multiqc report

# symbolically link clean reads to fastp_multiQC dir
# ln -s ../../../../../share/nefsc/mcfarland_sequecenes/TagSeq_oysters_2021/SA21200*/*.fastq.gz  ./ # call backward from the directory to the share folder, input symbolic link to folder
# hashed out if the symbolic links are already present

# Make an array of sequences to trim
array1=($(ls *.fastq.gz))  # call the folder will all symbolically linked .fastq.gz files (without the SA* folder included)

# fastqc loop of raw reads - output fastqc files to raw folder
for i in ${array1[@]}; do
        fastqc ${i} --outdir ./raw
done

echo "QC of raw reads complete." $(date)

# Quality Assessment of Raw reads

source ../../../python_venv/bin/activate # from the current directory, activates the bin of installed python packages, including multiqc

multiqc ./raw #Compile MultiQC report from FastQC files - output .html in current directory ( multiQC html report)

echo "Raw MultiQC report generated." $(date)


```

**Run the sbatch**

```
sbatch raw_multiQC
```

**Export multiqc report**

*exit bluewaves and run from terminal*
- save to gitrepo as multiqc_clean.html
```
scp samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/Geoduck_TagSeq/output/fastp_mutliQC/multiqc_report.html  C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/HPC_work/Output
```

### IMPORTANT! Quality check multiqc.html before you proceed!

- view the multiqc.html report and note observations to guide trimming!
  - **Per Sequence GC Content**: double peak of ~0-10% GC and 30-40% GC
**To do:** *initial peak due to polyA tail?*

  - **Mean Quality Scores*&: >22-2; **To do:** *trim base calls < 30*

  - **Per Base Sequence Content**: 15-25% C and G, 25-28% T, 35-40+% A

  - **Adapter Content**: high adapters present, 1-6% of sequences; **To do:** *trim adapter sequence*






# Trimming and post-trim quality check of 'clean' reads
-------------------------------------------











## fastp

**About:** preprocessing for FastQ files - open source code and info found on github [here](https://github.com/OpenGene/fastp)


### Trimming-polyA-tail
- Remember that TagSeq involves priming from the polyA tail of mRNA sequences! Thus, we will need to trim mononnucleotide sequence of As using fastp in addition to threshold quality score and adapters trimming. Furthermore, the default low complexity filter is 30% - given the prevalence of mononucleitide bases (i.e. polyA tail) it can be beneficial to increase this threshold to 50%


### What this script will do...
- ``` --adapter_sequence ``` =
	- trim adapter sequence ```AGATCGGAAGAGCACACGTCTGAACTCCAGTCA```
	- common single-end adapter in Illumina. You can run a test on a fastq.gz to count
- ``` --adapter_fasta``` = polyA_tail.fasta
	- **create 'polyA_tail.fasta' to call here**
	- important! ``` --adapter_sequence ``` is called by fastp before ``` --adapter_fasta``` and will call each adapter in the .fasta one-by-one
	- the sequence distribution of trimmed adapters can be found in the HTML/JSON report
- ```multiqc ./``` = outputs mutliqc report of the 'clean' reads in the current directory


### <span style="color:red">*sanity check!*

before we trim, trim, trim away, <span style="color:red">important to use stepwise diagnostics...

(1) first trim ONLY adapters and read the multiQC report

*if further processing is needed...*

(2) trim adapters + other parameters to optimize quality


# 'Adapter only' trim

# shell script: <span style="color:green">**fastp_multiQC_adapters_only**<span>

```
#!/bin/bash
#SBATCH --job-name="fastp_multiQC_adapters_only"
#SBATCH -t 002:00:00
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output=./Cvirginica_multistressor_TagSeq/output/fastp_multiQC/adapter_trim/"%x_out.%j"
#SBATCH --error=./Cvirginica_multistressor_TagSeq/output/fastp_multiQC/adapter_trim/"%x_err.%j"

# before running..
# make directory named adapter_trim in the folder output/fastp_multiQC/
# run sbatch from the host home directory

# load modules needed
module load bio/fastp/0.23.2
module load bio/fastqc/0.11.9

cd Cvirginica_multistressor_TagSeq/output/fastp_multiQC # run from the sgurr directory, navigate to the fast_multiQDC folder

# symbolically link clean reads to fastp_multiQC dir
# ln -s ../../../../../share/nefsc/mcfarland_sequecenes/TagSeq_oysters_2021/SA21200*/*.fastq.gz  ./ # call backward from the directory to the share folder, input symbolic link to folder
# commented out if the symbolic links are already created

# Make an array of sequences to trim
array1=($(ls *.fastq.gz))  # call the folder will all symbolically linked .fastq.gz files (without the SA* folder included)

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA)
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 ./adapter_trim/adapter_trim.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA # --trim_poly_x 6 -q 30 -y -Y 50  # check JUST adapters trimmed without the polyA tail and 50% complexity filter (-Y, defaults 30%)
        fastqc  ./adapter_trim/adapter_trim.${i} --outdir ./adapter_trim # call the output files from fastp in previous line and output fastqc in the same folder with adapter_trim filename head
done

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

source ../../../python_venv/bin/activate # from the current directory, activates the bin of installed python packages, including multiqc

multiqc ./adapter_trim  -o ./adapter_trim #Compile MultiQC report from FastQC files - output .html in adpater_trim directory ( fast_muiltiQC folder)

echo "Cleaned MultiQC report generated." $(date)

```

#### EXPORT MUTLIQC REPORT

*exit sedna and run from terminal*

- save to gitrepo  and rename as multiqc_report_adapter_trim_only.html

```
scp samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/Geoduck_TagSeq/output/clean/multiqc_report.html  C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics

```


### <span style="color:red">*report indicated further processing WAS needed...*



# 'Clean' trim (adapters + others)
# shell script: <span style="color:green">**fastp_multiQC_adapters_only**<span>

```
#!/bin/bash
#SBATCH --job-name="fastp_multiQC_clean"
#SBATCH -t 002:00:00
#SBATCH --mem=32GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output=./Cvirginica_multistressor_TagSeq/output/fastp_multiQC/clean/"%x_out.%j"
#SBATCH --error=./Cvirginica_multistressor_TagSeq/output/fastp_multiQC/clean/"%x_err.%j"

# before running..
# make directory named clean in the folder output/fastp_multiQC/
# run sbatch from the host home directory

# load modules needed
module load bio/fastp/0.23.2
module load bio/fastqc/0.11.9

cd Cvirginica_multistressor_TagSeq/output/fastp_multiQC # run from the sgurr directory, navigate to the fast_multiQDC folder

# symbolically link clean reads to fastp_multiQC dir
# commented out if the symbolic links are already created
# ln -s ../../../../../share/nefsc/mcfarland_sequecenes/TagSeq_oysters_2021/SA21200*/*.fastq.gz  ./ # call backward from the directory to the share folder, input symbolic link to folder

# Make an array of sequences to trim
array1=($(ls *.fastq.gz))  # call the folder will all symbolically linked .fastq.gz files (without the SA* folder included)

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA)
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 ./clean/clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50
        fastqc  ./clean/clean.${i} --outdir ./clean # calls the  files output by fastp in previous line and output into the same folder
done

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

source ../../../python_venv/bin/activate # from the current directory, activates the bin of installed python packages, including multiqc

multiqc ./clean -o ./clean #Compile MultiQC report from FastQC files - output .html in current directory ( fast_muiltiQC folder)

echo "Cleaned MultiQC report generated." $(date)

```

### EXPORT MUTLIQC REPORT

*exit sedna and run from terminal*

- save to gitrepo and rename as multiqc_report_clean.html

```
scp samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/Geoduck_TagSeq/output/clean/multiqc_report.html  C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics

```


## we now have **'clean' fastq files** after processing with fastp! Use these in the following...









# Alignment of cleaned reads to reference
-------------------------------------------








## Reference genome upload to HPC

*exit HPC and run from terminal*

-	move refs Crassostrea virginica to the HPC, references include the following (downloaded from https://www.ncbi.nlm.nih.gov/genome/?term=txid6565[Organism:noexp]):

* genome (genomic.fna)
* transcript reference (rna.fna)
* anotation (.gff)   

### import species refs to the HPC

Below is the scp calling each reference .gz file based on common GCF filename header - import to the HPC

```
scp C:/Users/samuel.gurr/Documents/Genomes/Crassostrea_virginica/GCF* sgurr@sedna.nwfsc2.noaa.gov:refs/
```

### gunzip within the HPC

*ssh back into the HPC and navigate to the folder with uplaoded refs*

```
gunzip GCF*
```

now 'ls' and see that .gz is gone, all files are unzipped and ready for next steps!


### citations for *C. virginica* references

-  file name: GCF_002022765.2_C_virginica-3.0_genomic.fna
	- Gómez-Chiarri M et al., "Developing tools for the study of molluscan immunity: The sequencing of the genome of the eastern oyster, Crassostrea virginica.", Fish Shellfish Immunol, 2015 Sep;46(1):2-4
	- Milbury CA et al., "Complete mitochondrial DNA sequence of the eastern oyster Crassostrea virginica.", Mar Biotechnol (NY), 2005 Nov-Dec;7(6):697-712
-  reference genome file size: total length (Mb) 684.741

### **now we are ready!!**


## HISAT2 alignment

**About:** HISAT2 is a sensitive alignment program for mapping sequencing reads.
In our case here, we will use HISAT2 to **(1)** index our *P. generosa* reference genome **(2)** align our clean TagSeq reads to the indexed reference genome.

More information on HISAT2 can be read [here](http://daehwankimlab.github.io/hisat2/manual/)!

*Main arguments used below*:

**(1)** *Index the reference genome*

``` hisat2-build ``` =
builds a HISAT2 index from a set of DNA sequences. Outputs 6 files that together consitute the index.
ALL output files are needed to asslign reads to the reference genome and the original sequence FASTA file(s)
are no longer used for th HISAT2 alignment

``` -f <reads.fasta> ``` =
the reads (i.e <m1>, <m2>, <m100>)
FASTA files usually have extension .fa, .fasta, .mfa, .fna or similar.
FASTA files do not have a way of specifying quality values, so when -f is set,
the result is as if --ignore-quals is also set

Note: other options for your reads are ```-q ``` for FASTQ files,
 ```--qseq ``` for QSEQ  files, etc. - check [here](http://daehwankimlab.github.io/hisat2/manual/) for more file types

**(2)** *align reads to reference*

``` -x <hisat2-indx> ``` =
base name for the index reference genome, first looks in the current directory then in the directory specified in HISAT_INDEXES environment variable

``` --dta ``` =
 **important!** reports the alignments tailored for *StringTie* assembler.
 With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites.
 This leads to fewer alignments with short-anchors, which helps transcript assemblers improve
 significantly in computation and memory usage.

 This is important relative to Cufflinks ```--dta-cufflinks```
 in which HISAT2 looks for novel splice sites with three signals (GT/AG, GC/AG, AT/AC)

``` -U <r> ``` =
Comma-separated list of files contained unparied reads to be aligned.
*In other words...*, our array of post-trimmed 'clean' TagSeq reads

``` -p NTHREADS``` =
Runs on separate processors, increasing -p increases HISAT2's memory footprint, increasing -p from 1 to 8
increased the footprint by a few hundred megabytes

Note: Erin Chile from Putnam Lab ran -p 8 for HISAT2 alignment of *Montipora* [here](https://github.com/echille/Montipora_OA_Development_Timeseries/blob/master/Amil/amil_RNAseq-analysis.sh)

``` -S <hit> ``` =
file to write SAM alignments to. By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console).

**Note:** HISAT2 also has several criteria to trim such as the phred score and base count 5' (left) and 3' (right)
Since we already assessed quality and trimmed, we will not use these commands


## samtools
**About:** used to manipuate alignments to a binary BAM format - files contain the spliced reads alignemnts sorted by the reference position
with a tag to indicate the genomic strand that produced the RNA from which the read was sequenced.
samtools quickly extracts alignments overlapping particular genomic regions - outputs allow viewers to quickly display alignments in each genomic region
note: the SAM format output by HISAT2 must be sorted and converted to BAM format using the samtools program

``` sort ``` =
sorts the alignments by the leftmost coordinates

```-o <out.bam``` =
outputs as a specified .bam file

```-@ threads``` =
Set number of sorting and compression threads. By default, operation is single-threaded

more information on samtools commands [here](http://www.htslib.org/doc/1.1/samtools.html)


## python3
hisat2-build has contingencies for python
python3 is located in the virtual envrionment so we will need to call an alias for it in our .bashrc
hisat-build should now access python in the following script

```
cd

nano .bashrc

alias python=/python_env/bin/python3

```



## HPC Job: HISAT2 Index Reference and Alignment
-----------------------------------------------------------------

- index reference and alignment

- create directory output\hisat2

``` mkdir hisat2 ```

**input**
- GCF_002022765.2_C_virginica-3.0_genomic.fna *= reference genome*
- clean/*.fastq.gz *= all clean TagSeq reads from fastp*

**ouput**
- Cvirginica_ref *= indexed reference by hisat2-build; stored in the output/hisat2 folder as 1.hy2, 2.ht2... 8.ht2*
- <clean.fasta>.sam *=hisat2 output, readable text file; removed at the end of the script*
- <clean.fasta>.bam *=converted binary file complementary to the hisat sam files*

# shell script: <span style="color:green">**hisat2.sh**<span>

```
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


```
### HISAT2 complete with format **prepared for StringTie assembler!**


### <span style="color:red">IMPORTANT:<span> merge lanes here (.bam stage)

NOTE: you may see in the seq filenames something  like L001 and L002 indicating that the same sample was sequenced on two lanes. You can proceed here to combine the bam files together using samtools - this does not interfere with the gene count matrix (I have tested both merging and without merging)

Why merge? - the count matrix will have column for each filename, the merge simply reduces a single step in R to sum together columns with same sample ID on multiple lanes.

If you want to proceed with the merge....

- enter interactive mode

```
interactive
```

- load samtools
```
module load SAMtools/1.9-foss-2018b
```

- ID file lists the sample ID characters (i.e. SG9_S108) - sample .bam file SG98_S144_L002_R1_001.bam
```
ls *R1_001.bam | awk -F '[_]' '{print $1"_"$2}' | sort | uniq > ID
```

- for loop using ```samtools``` to merge bam files together by sample ID ('ID' created above)
- example of files to merge: SG98_S144_L001_R1_001.bam and SG98_S144_L002_R1_001.bam - where 'L001' and 'L002' are the lanes
and the ID called SG9_S108 as $i below... (note: all .gz file names end with '001.fasta.gz' no matter the lane)
```
for i in `cat ./ID`;
	do samtools merge $i\.bam $i\_L001_R1_001.bam $i\_L002_R1_001.bam;
	done
```

- when run in ```interactive``` mode, this loop takes up to ~1-2 hours

# Assembly and quantification
-----------------------------------------------------------------

### Upload annotation reference .gff or .gff3 to HPC

- Sam White and Steven Roberts completed the annotation with open resources available for download

*exit bluewaves and run from terminal*

- "one file to rule them all" (- Sam W.) with the CDS, mRNA, exon, and genes
```
scp C:/Users/samjg/Documents/Bioinformatics/genomes/Pgenerosa_annotations/Panopea-generosa-vv0.74.a3-merged-2019-09-03-6-14-33.gff3 samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/refs/
```

- mRNA only
```
scp C:/Users/samjg/Documents/Bioinformatics/genomes/Pgenerosa_annotations/Panopea-generosa-v1.0.a4.gene.gff3 samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/refs/
```

- mRNA annotation with 'm01, m02, m03, m04... m10' removed from gene.ID to match IDs

```
scp C:/Users/samjg/Documents/Bioinformatics/genomes/Pgenerosa_annotations/Panopea-generosa-v1.0.a4.mRNA_SJG.gff3 samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/refs/
```

## StringTie
-----------------------------------------

**About:** StringTie is a fast and efficient assembler of RNA-Seq alignments to assemble and quantitate full-length transcripts.  putative transcripts.
For our use, we will input short mapped reads from HISAT2. StringTie's ouput can be used to identify DEGs in programs such as DESeq2 and edgeR

More information on StringTie can be read [here](https://ccb.jhu.edu/software/stringtie/)!

[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) is another assembler option we will not use. [Read this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4643835/#:~:text=On%20a%20simulated%20data%20set,other%20assembly%20software%2C%20including%20Cufflinks.)

*Main StringTie arguments used below*:

``` -p ``` =
specify the number of threads (CPUs). Note that a single node has at least 8 CPUs or 'threads' to call here

``` -A ``` =
Gene abundances will be reported (tab delimited format) in the output file with the given name.

``` -e ``` =
Limits the processing of read alignments to only estimate and output the assembled transcripts matching the
reference transcripts given with the -G option (*requires -G, recommended for -B/-b*). **With this option,
read bundles with no reference transcripts will be entirely skipped**, which may provide a considerable
speed boost when the given set of reference transcripts is limited to a set of target genes, for example.

``` -G <ref_ann.gff>``` =
Use the reference annotation file (in GTF or GFF3 format) to guide the assembly process.
The output will include expressed reference transcripts as well as any novel transcripts that are assembled.
*This option is required by options -B, -b, -e, -C*.

``` -o ``` =
Sets the name of the output GTF file where StringTie will write the assembled transcripts

``` --merge ``` =
StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. This mode is used in the new differential analysis pipeline to generate a global, unified set of transcripts (isoforms) across multiple RNA-Seq samples.
If the -G option (reference annotation) is provided, StringTie will assemble the transfrags from the input GTF files with the reference transcripts.

- merge commands

 ```-G <guide_gff>```	reference annotation to include in the merging (GTF/GFF3)

 ```-o <out_gtf>```	output file name for the merged transcripts GTF (default: stdout)

## gffcompare
-----------------------------

**About:** gffcompare is used to estimate the accuracy of one or more GFF query files compared to a reference annotation

*Main gffcompare arguments used below*:

``` -r ``` =
reference annotation GFF file - each file is matched against it and tagged as overlapping, matching, or novel where appropriate.
outputs .refmap and .tmap files

``` -G ``` =

``` -o <outprefix>``` = give a prefix for output files created by gffcompare (i.e. merged)


## Python step prepDE.py
-----------------------------

**About:** in preparation for differential expression analysis using Bioconductor packages in R (i.e. DESeq2, edgeR, WGCNA),
we will need to first aquire a matrix of read counts to particular genomic features (i.e. genes). prepDE.py is a Python call to extract
hypothetical read counts for each transcript from the files generated by ```-e``` using StringTie (here as ${i}.gtf from .bam files by HISAT2)

```prepDE.py``` =
.py found [here](https://github.com/gpertea/stringtie/blob/master/prepDE.py) - add to scripts folder to call using python
prepDE.py, builds the matrix of read counts

``` -g ``` =
where to ouput the gene count matrix (defaults as gene_count_matrix.csv)
```-i INPUT``` =
a folder containing all sample sub-directories; Alternatively can provide a text file (i.e. sample_list.tct) with sample ID and path to its
GTF file on each line (default '.' meaning the working directory is the subdirectory with all GTF files)
Alternatively call a **listGTF.txt** file... this file has two columns with the sample ID and the <Path to sample> for the .gtf files
run the following:

* move prepDE.py (downloaded onlne) to the ssh 

```
scp C:/Users/samuel.gurr/Documents/Github_repositories/Cvriginica_multistressor/HPC_Analysis/Scripts/ prepDE.py sgurr@sedna.nwfsc2.noaa.gov:Cvirginica_multistressor_TagSeq/scripts/
```
-----------------------------

**gtf_list.txt** run...

- ``` ls .gtf > gtf_list.txt```
- ```--merge``` requires a list file to call each of the gtf files

**listGTF.txt** run...

- ```for filename in *.gtf; do echo $filename $PWD/$filename; done > listGTF.txt```
- call this text file ``` -i ``` in  python prepDE.py in the job merge_prepDEpy.sh


# HPC job: Assembly
-----------------------------------------------------------------

# shell script: <span style="color:green">**stringtie2.sh**<span>
```
#!/bin/bash
#SBATCH --job-name="stringtie2"
#SBATCH -t 048:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"

# load modules, requires stringtie2
module load bio/stringtie/2.2.0

# before running..
# make directory named stringtie2 as output/stringtie2/

cd #nav back to home directory (allows job to be run from anywhere)
cd Cvirginica_multistressor_TagSeq/output/stringtie2 #nav back to target directory stringtie2

array=($(ls ../hisat2/merged/*.bam)) #Make an array of sequences to assemble

# awk example: file name clean.C9-larva-22_S77.bam - use awk to get sample name w/o .bam using . as delimiter 
# run strnigtie on the array and output to the stringtie2 directory
for i in ${array[@]}; do #Running with the -e option to compare output to exclude novel genes. Also output a file with the gene abundances
        sample_name=`echo $i| awk -F [.] '{print $1"_"$2}'`
	stringtie -p 8 -e -B -G ../../refs/GCF_002022765.2_C_virginica-3.0_genomic.gff -A ../../stringtie2/${sample_name}.gene_abund.tab -o ../../stringtie2/${sample_name}.gtf ${i}
        echo "StringTie assembly for seq file ${i}" $(date)
done
echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)

```





# HPC job: Merge and Build Read Count Matrix for DEG analysis
-----------------------------------------------------------------

- NOTE: you will need the files **gtf_list.txt** and **listGTF.txt** to in your -D working directory (i.e. output/stringtie) to run this job (described above)

# shell script: <span style="color:green">**stringtie2_merge_prepDEpy.sh**<span>
```
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
#SBATCH -D /data/putnamlab/sgurr/Geoduck_TagSeq/output/stringtie

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
```
*exit bluewaves and run from terminal*

- save the read count matrix to local pc
```
scp  samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/Geoduck_TagSeq/output/stringtie/lanes_merged/transcript_count_matrix.csv C:/Users/samjg/Documents/My_Projects/Pgenerosa_TagSeq_Metabolomics/TagSeq/HPC_Bioinf
```
