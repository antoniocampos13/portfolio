# FastQ to Annotation bash script

Variant annotation pipeline intended for automated next generation DNA sequencing analysis. It takes next generation sequencing (NGS) raw reads in the form of gzip-compressed FASTQ files (.fast.gz) from human whole genome sequencing as input and produces variant annotation as output. Variant annotation is the process of identifying genetic variants in some genomic DNA sample, and assess, for example, if any of the found variants have any effect on phenotype, such as increased susceptibility to certain diseases.

This script is a simplified, modified version of the work by Blighe, Beauchamp, and colleagues at Sheffield Diagnostic Genetics Service, Sheffield Children's NHS Foundation Trust, Sheffield, UK, and their efforts to introduce a clinical-grade next generation sequencing (NGS) analysis pipeline fully validated against Sanger di-deoxy sequencing.

The pipeline is built using open source programs mixed with customised scripts. A master and concise log is kept, with date- and time-stamps. Results directory structure is formed based on the filename of mate-pair `FASTQ` files.

The unique feature of the analysis pipeline that increases sensitivity to Sanger sequencing is in the variant calling step, where a final aligned BAM is split into 3 'sub-BAMs' of 75%, 50%, and 25% random reads. Variants are then called on all 4 BAMs, after which a consensus VCF is produced.

This pipeline proceeds in an 8-step process:

* Adaptor and read quality trimming - TrimGalore! (Krueger F), FastQC (Andrews S), cutadapt (Martin M, 2011)
* Alignment - bwa mem (Li & Durbin, 2009)
* Marking and removing PCR duplicates - Picard (Broad Institute of MIT and Harvard), SAMtools (Li et al., 2009)
* Remove low mapping quality reads - SAMtools (Li et al., 2009)
* QC - SAMtools (Li et al., 2009), BEDTools (Quinlan & Hall, 2010), custom scripts
* Downsampling/random read sampling - Picard (Broad Institute of MIT and Harvard)
* Variant calling - SAMtools/BCFtools (Li et al., 2009)
* Annotation - Variant Effect Predictor (McLaren et al., 2016)

## Requirements

* (mini)conda environment with necessary programs. Recommended installation:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -y --name <YOUR_ENV> python=3.6
conda activate <YOUR_ENV>
conda install --yes --file requirements.txt
```

* Local RefSeq cache of Ensembl Variant Effect Predictor (VEP):

```bash
vep_install -a cf -s homo_sapiens_refseq -y GRCh38 -c . –CONVERT

# Alternative:
wget ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_refseq_vep_101_GRCh38.tar.gz -P $HOME/.vep

tar -zxf $HOME/.vep/homo_sapiens_refseq_vep_101_GRCh38.tar.gz -C $HOME/.vep
```

* Paired-end gz-compressed `FASTQ` files
* Human genome reference `FASTA` file and corresponding `bwa` indexes, as well as chromosome-ordered and position-sorted BED or GFF file in GRCh38 format. Recommended setup:

```bash
# Download GRCh38 major release without ALT contigs and with decoy genomes (EBV and hs38d1 contig) from NCBI's FTP server
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

# Uncompress file and then perform bwa indexing
bwa index GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

# Produce samtools index
samtools faidx GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

# Or download from NCBI's FTP server
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai

# Produce sorted BED file from reference genome index file
awk '{print $1 "\t0\t" $2}' GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai | sort -k1,1V -k2,2n > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.bed

# Or use RefSeq GFF file instead of sorted BED file
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
```

## Execution

Activate your conda environment, then run the script using the following parameters:

* Mate-pair FASTQ files name root (without extension) (absolute file path)
* Reference genome FASTA (absolute file path)
* BED or GFF file (absolute file path)
* Minimum quality for bases at read ends, below which bases will be cut (integer - default: 20)
* Minimum allowed read length (integer - default: 20)
* Adaptor for trimming off read ends ('illumina' / 'nextera' / 'small_rna')
* Minimum read depth for calling a variant (integer - default:)
* Minimum allowed mapping quality (integer - default: )
* Stringency for calling variants ('relaxed' / 'normal') (relaxed uses --pval-threshold 1.0 with BCFtools call)
* User identification for logging (alphanumeric)

The script include a preliminary step that will check if all parameters are adequate.

## Output

Output files are generated at the directory where input FASTQ files are. The final part of the script is a tidy-up step that organizes all generated files to a directory named  prefix_results where the prefix is the FASTQ mate-pair accession ID. This folder will contain:

* Log files
  * prefix_Pipeline_Log.txt - short log
  * prefix_Master_Log.txt - comprehensive log
  
* A variant call file 
  * prefix_Final.vcf

* Annotation files
  * prefix_AnnotationVEP.html - HTML report of variant annotation, with consequences for all known transcript isoforms
  * prefix_AnnotationVEP.txt - as above but in tab-separated text format

Additionally, other three subfolders will be created to organize intermediate files:

* trimmed_files folder
  * Trimmed read file (fastq.gz) for each mate-pair
  * prefix_trimming_report.txt for each mate-pair
  * FastQC reports (in HTML format) for each mate-pair

* alignment_files folder
  * prefix_Alignment.txt - alignment metrics
  * prefix_ReadsOffTarget.txt - number of reads falling outside regions specified in BED/GFF file
  * prefix_PCRDuplicates.txt - details on identified PCR duplicates
  * prefix_CoverageTotal.bedgraph - coverage for all mapped locations (contiguous bases at same read depth are merged into regions)
  * prefix_MeanCoverageBED.bedgraph - mean read depth for each region specified in supplied BED/GFF file
  * prefix_PerBaseDepthBED.bedgraph - per base read depth for each base in each region specified in supplied BED/GFF file
  * prefix_PercentGenomeCovered.txt - percentage of reference genome covered by reads.
  * prefix_Aligned_Sorted_PCRDuped_FiltMAPQ.bam - aligned BAM file with sorted reads, PCR duplicates removed, and reads below mapping quality threshold removed
  * prefix_Aligned_Sorted_PCRDuped_FiltMAPQ.bam.bai - index for above BAM file

* variant_call_files folder
  * samples.txt - file created to assist downsampling/random read sampling step
  * Intermediate VCF files

### Original Authors

* Kevin Blighe (Sheffield Children's NHS Foundation Trust)
* Nick Beauchamp (Sheffield Children's NHS Foundation Trust)
* Darren Grafham (Sheffield Children's NHS Foundation Trust)
* Lucy Crookes (Sheffield Children's NHS Foundation Trust)
* Sasirekha Palaniswamy ('Sashi') (Sheffield Children's NHS Foundation Trust)
* Sheffield Diagnostics Genetics Service

### Author of the present (modified) version

* Antonio Victor Campos Coelho ([Contact me](https://antoniocampos13.github.io/pages/contact.html#contact))

# References

[Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Picard Tools - By Broad Institute](http://broadinstitute.github.io/picard/)

[Babraham Bioinformatics - Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Li H and Durbin R (2009), Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics 25(14): 1754–1760.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup (2009), The Sequence Alignment/Map format and SAMtools, Bioinformatics 25(16): 2078-9.

Martin M (2011), Cutadapt removes adapter sequences from high-throughput sequencing reads, EMBnet.journal 17(1): 10-12.

McLaren W, Gil L, Hunt S, Riat H, Ritchie G, Thormann A, Flicek P, Cunningham F (2016), The Ensembl Variant Effect Predictor, Genome Biology 17: 122.

Quinlan AR & Hall IM (2010), BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics 26(6): 841-2.