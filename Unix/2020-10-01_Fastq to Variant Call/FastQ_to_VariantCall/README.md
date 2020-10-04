# FastQ to Variant Call bash script

Variant annotation pipeline intended for automated next generation DNA sequencing analysis. It takes next generation sequencing (NGS) raw reads in the form of gzip-compressed FASTQ files (.fast.gz) from human whole genome sequencing as input and produces variant annotation as output. Variant annotation is the process of identifying genetic variants in some genomic DNA sample, and assess, for example, if any of the found variants have any effect on phenotype, such as increased susceptibility to certain diseases.

This script is a simplified, modified version of the work by Blighe, Beauchamp, and colleagues at Sheffield Diagnostic Genetics Service, Sheffield Children's NHS Foundation Trust, Sheffield, UK, and their efforts to introduce a clinical-grade next generation sequencing (NGS) analysis pipeline fully validated against Sanger di-deoxy sequencing.

The pipeline is built using open source programs mixed with customised scripts. A master and concise log is kept, with date- and time-stamps. Results directory structure is formed based on the run number and patient ID.

The unique feature of the analysis pipeline that increases sensitivity to Sanger sequencing is in the variant calling step, where a final aligned BAM is split into 3 'sub-BAMs' of 75%, 50%, and 25% random reads. Variants are then called on all 4 BAMs, after which a consensus VCF is produced.

This pipeline proceeds in an 8-step process:

* Adaptor and read quality trimming - TrimGalore! (Krueger F), FastQC (Andrews S), cutadapt (Martin M, 2011)
* Alignment - bwa mem (Li & Durbin, 2009)
* Marking and removing PCR duplicates - Picard (Broad Institute of MIT and Harvard), SAMtools (Li et al., 2009)
* Remove low mapping quality reads - SAMtools (Li et al., 2009)
* QC - SAMtools (Li et al., 2009), BEDTools (Quinlan & Hall, 2010), custom scripts
* Downsampling / random read sampling - Picard (Broad Institute of MIT and Harvard)
* Variant calling - SAMtools/BCFtools (Li et al., 2009)
* Annotation - Variant Effect Predictor (McLaren et al., 2016)

## Requirements

* (mini)conda environment with necessary programs
  * Create fastqvc environment with miniconda_requirements.yml
  * You may change the environment name in the yml file if desired

* (Optional, but recommended) Local RefSeq cache of Ensembl Variant Effect Predictor (VEP) installed via conda environment:
  * vep_install -a cf -s homo_sapiens_refseq -y GRCh38 -c . –CONVERT

* Paired-end gz-compressed FASTQ files
* Chromosome-ordered and position-sorted BED or GFF file in hg38

## Execution

Run the script, which will check command-line parameters, execute the main pipeline, and then return results files. Use the following parameters:

* FASTQ mate-pair accession ID (absolute file path)
* Reference genome FASTA (absolute file path)
* Run number (e.g. Plate6, Plate7, etc.) (alphanumeric)
* Patient ID (alphanumeric)
* BED or GFF file (absolute file path)
* Minimum quality for bases at read ends, below which bases will be cut (integer - default: 20)
* Minimum allowed read length (integer - default: 20)
* Adaptor for trimming off read ends ('illumina' / 'nextera' / 'small_rna')
* Minimum read depth for calling a variant (integer)
* Minimum allowed mapping quality (integer)
* Stringency for calling variants ('relaxed' / 'normal') (relaxed uses --pval-threshold 1.0 with BCFtools call)
* Directory where results will be output (absolute file path)
* User initials (alphanumeric)

### Hard-coded sections of code

* FastQ_to_VariantCall.sh, line 32: fastqvc miniconda environment. Change as needed.

## Output

Results files are output locally to [results root]/[run number]/[sample ID]/
*_AnalysisLog.txt - analysis log (short)
* Master.log - analysis log (comprehensive)
*_R1_001.fastq.gz_trimming_report.txt - details on base and read trimming for mate-pair 1
*_R1_001_val_1_fastqc.html - FastQC report for mate-pair 1 (after trimming)
*_R2_001.fastq.gz_trimming_report.txt - details on base and read trimming for mate-pair 2
*_R2_001_val_2_fastqc.html - FastQC report for mate-pair 2 (after trimming)
*_Alignment.txt - alignment metrics
*_ReadsOffTarget.txt - number of reads falling outside regions specified in BED file
*_PCRDuplicates.txt - details on identified PCR duplicates
*_CoverageTotal.bedgraph - coverage for all mapped locations (contiguous bases at same read depth are merged into regions)
*_MeanCoverageBED.bedgraph - mean read depth for each region specified in supplied BED file
*_PerBaseDepthBED.bedgraph - per base read depth for each base in each region specified in supplied BED file
*_PercentGenomeCovered.txt - percentage of reference genome covered by reads.
*_Aligned_Sorted_PCRDuped_FiltMAPQ.bam - aligned BAM file with sorted reads, PCR duplicates removed, and reads below mapping quality threshold removed
*_Aligned_Sorted_PCRDuped_FiltMAPQ.bam.bai - index for above BAM file
*_Final.vcf - final VCF file
*_AnnotationVEP.html - HTML report of variant annotation, with consequences for all known transcript isoforms
*_AnnotationVEP.tsv - as above but in tab-separated values (TSV) format

### Original Authors

* Kevin Blighe (Sheffield Children's NHS Foundation Trust)
* Nick Beauchamp (Sheffield Children's NHS Foundation Trust)
* Darren Grafham (Sheffield Children's NHS Foundation Trust)
* Lucy Crookes (Sheffield Children's NHS Foundation Trust)
* Sasirekha Palaniswamy ('Sashi') (Sheffield Children's NHS Foundation Trust)
* Sheffield Diagnostics Genetics Service

### Author of the present (modified) version

* Antonio Victor Campos Coelho ([Contact me](https://antoniocampos13.github.io/pages/contact.html#contact))

# Reference
References

[Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Picard Tools - By Broad Institute](http://broadinstitute.github.io/picard/)

[Babraham Bioinformatics - Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Li H and Durbin R (2009), Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics 25(14): 1754–1760.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup (2009), The Sequence Alignment/Map format and SAMtools, Bioinformatics 25(16): 2078-9.

Martin M (2011), Cutadapt removes adapter sequences from high-throughput sequencing reads, EMBnet.journal 17(1): 10-12.

McLaren W, Gil L, Hunt S, Riat H, Ritchie G, Thormann A, Flicek P, Cunningham F (2016), The Ensembl Variant Effect Predictor, Genome Biology 17: 122.

Quinlan AR & Hall IM (2010), BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics 26(6): 841-2.