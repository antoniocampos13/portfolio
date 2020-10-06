#!/usr/bin/env bash

# Part 1 - Searching NCBI's databases for relevant projects
## Activate the environment
conda activate bioenv

## Let's create a folder to organize our files inside
mkdir demo
cd demo

## Basic search
esearch -db bioproject -query 'vorinostat' | efetch -format native -mode xml | xtract -pattern DocumentSummary -element ArchiveID@accession ID Reference Title Description > vorinostat_projects.txt

## Refined search
esearch -db bioproject -query '(((vorinostat) AND "bioproject sra"[Filter]) AND Homo sapiens[Organism]) AND "variation"[Filter]'  | efetch -format native -mode xml | xtract -pattern DocumentSummary -element ArchiveID@accession ID Reference Title Description  > vorinostat_refined.txt

# Part 2 - Connecting to SRA to obtain FASTQ files
## Conect to SRA database and retrieve information with a BioProject accession ID
esearch -db sra -query 'PRJNA436005' | efetch -format runinfo > PRJNA436005_runinfo.csv

## Download FASTQ files from project deposited in BioProject database using read set accession ID
fastq-dump --split-files SRR6784104 --gzip

# Part 3 - Final preparations
## Download and install VEP cache to $HOME/.vep
vep_install -a cf -s homo_sapiens_refseq -y GRCh38 -c . –CONVERT

### Alternative: manually download and uncompress the cache
wget ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_refseq_vep_101_GRCh38.tar.gz -P $HOME/.vep

tar -zxf $HOME/.vep/homo_sapiens_refseq_vep_101_GRCh38.tar.gz -C $HOME/.vep

## Download/generate human genome reference and index files
### First, create a subfolder into a demo folder to better organize our reference files
mkdir refs
cd refs

### 1. Human genome FASTA: Download GRCh38 major release without ALT contigs and with decoy genomes (EBV and hs38d1 contig) from NCBI's FTP server
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

### 2. Human genome FASTA index
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai

### 3. Human genome BED file: produce sorted BED file from reference genome index file obtained above
awk '{print $1 "\t0\t" $2}' GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai | sort -k1,1V -k2,2n > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.bed

### 4. Human genome GFF file (optional alternative to BED file)
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz

#### 5. bwa index files
# Uncompress file obtained with step 1
bwa index GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

# Part 4 - Running the script
cd .. # go back to demo folder

REF=refs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
BED=refs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bed

./FastQ_to_Annotation.sh SRR6784104 $REF $BED 20 20 illumina 3 0 normal your_name
