#!/usr/bin/env bash

#########################################################
#Pre-liminary step:  Safeguarding the script's operation#
#########################################################

echo "Checking command line parameters..."

Read1="/mnt/c/Users/victo/Documents/ubuntu/ClinicalGradeDNAseq/reads/${1}_1.fastq.gz"
Read2="/mnt/c/Users/victo/Documents/ubuntu/ClinicalGradeDNAseq/reads/${1}_2.fastq.gz"
Ref_FASTA="/mnt/c/Users/victo/Documents/ubuntu/ClinicalGradeDNAseq/refs/${2}"
RunNumber="${3}"
PatientID="${4}"
BEDorGFFfile="/mnt/c/Users/victo/Documents/ubuntu/ClinicalGradeDNAseq/refs/${5}"
TrimmingQualityReadEnds="${6}"
TrimmingReadLengthMin="${7}"
TrimmingAdaptor="${8}"
FilterReadDepthCutoff="${9}"
FilterMappingQualityCutOff="${10}"
CallingStringency="${11}"
Results_root="/mnt/c/Users/victo/Documents/ubuntu/ClinicalGradeDNAseq/${12}"
User="${13}"

#Check that 13 parameters have been passed
if [ $# -ne 13 ]
then
        echo "Error - incorrect number of parameters. Consult the relevant standard operating procedure for correct usage."
	exit 1
fi

#Check that the file $Read1 exists
test -e "${Read1}"
if [ $? -ne 0 ]
then
        echo "Error - ${Read1} does not exist. Please check the complete file path and re-run."
	exit 1
fi

#Check that the file $Read2 exists
test -e "${Read2}"
if [ $? -ne 0 ]
then
        echo "Error - ${Read2} does not exist. Please check the complete file path and re-run."
        exit 1
fi

#Check that the file $Ref_FASTA exists
test -e "${Ref_FASTA}"
if [ $? -ne 0 ]
then
        echo "Error - ${Ref_FASTA} does not exist. Please check the complete file path and re-run."
        exit 1
fi

#Check that the file $BEDorGFFfile exists
test -e "${BEDorGFFfile}"
if [ $? -ne 0 ]
then
        echo "Error - ${BEDorGFFfile} does not exist. Please check the complete file path and re-run. GFF files are also accepted."
        exit 1
fi

#Check that various parameters are integers
if [[ ! "${TrimmingQualityReadEnds}" =~ ^[0-9]*$ ]]
then
	echo "Error - "${TrimmingQualityReadEnds}" must be an integer."
	exit 1
fi
if [[ ! $TrimmingReadLengthMin =~ ^[0-9]*$ ]]
then
        echo "Error - "${TrimmingReadLengthMin}" must be an integer."
        exit 1
fi
if [[ ! $FilterReadDepthCutoff =~ ^[0-9]*$ ]]
then
        echo "Error - "${FilterReadDepthCutoff}" must be an integer."
        exit 1
fi
if [[ ! $FilterMappingQualityCutOff =~ ^[0-9]*$ ]]
then
        echo "Error - "${FilterMappingQualityCutOff}" must be an integer."
        exit 1
fi

#Check that the specified trimming adaptor is correct
if [ $TrimmingAdaptor != "illumina" ]
then
	if [ $TrimmingAdaptor != "nextera" ]
	then
		if [ $TrimmingAdaptor != "small_rna" ]
		then
		        echo "Error - "${TrimmingAdaptor}" must be one of illumina|nextera|small_rna."
		        exit 1
		fi
	fi
fi

#Check that the specified trimming adaptor is correct
if [ $CallingStringency != "relaxed" ]
then
        if [ $CallingStringency != "normal" ]
        then
		echo "Error - "${CallingStringency}" must be one of normal|relaxed."
		exit 1
        fi
fi

echo "Done."

echo "Running master analysis script on `date`..."

mkdir -p "${Results_root}"
mkdir -p "${Results_root}"/"${RunNumber}"
mkdir -p "${Results_root}"/"${RunNumber}"/"${PatientID}"
mkdir -p "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp

echo "Log file is "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${4}"_Master.log"

eval "$(conda shell.bash hook)"
conda activate bioenv

# echo -e "\n\n\n"
# echo "###########################################################################"
# echo "#Pre-liminary step:  Log initialisation#"
# echo "###########################################################################"

#Begin log file
echo "Beginning analysis script on `date`, run by "${User}" with the following parameters:" \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t1\t"${Read1}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t2\t"${Read2}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t3\t"${Ref_FASTA}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t4\t"${RunNumber}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t5\t"${PatientID}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t6\t"${BEDorGFFfile}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t7\t"${TrimmingQualityReadEnds}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t8\t"${TrimmingReadLengthMin}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t9\t"${TrimmingAdaptor}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t10\t"${FilterReadDepthCutoff}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t11\t"${FilterMappingQualityCutOff}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t12\t"${CallingStringency}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t13\t"${Results_root}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt
echo -e "\t14\t"${User}"" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

echo "Done."

exec &> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${4}"_Master.log

echo -e "\n\n"
echo "#####################################################"
echo "#Analysis step 1:  adaptor and read quality trimming#"
echo "#####################################################"
echo "Beginning analysis step 1 (adaptor and read quality trimming) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

trim_galore \
	--qual "${TrimmingQualityReadEnds}" \
	--gzip \
	--length "${TrimmingReadLengthMin}" \
	--"${TrimmingAdaptor}" \
	--paired "${Read1}" "${Read2}" \
	--output_dir "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/ \
	--fastqc

TrimmedPath=""${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"
TrimmedFile1="${TrimmedPath}""${1}"_1_val_1.fq.gz
TrimmedFile2="${TrimmedPath}""${1}"_2_val_2.fq.gz
TrimmingReport1="${TrimmedPath}""${1}"_1.fastq.gz_trimming_report.txt
TrimmingReport2="${TrimmedPath}""${1}"_2.fastq.gz_trimming_report.txt
TrimmingHtml1="${TrimmedPath}""${1}"_1_val_1_fastqc.html
TrimmingHtml2="${TrimmedPath}""${1}"_2_val_2_fastqc.html

echo "Done."



echo -e "\n\n"
echo "#############################"
echo "#Analysis step 2:  Alignment#"
echo "#############################"
echo "Beginning analysis step 2 (alignment) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

#BWA mem
bwa mem \
	"${Ref_FASTA}" \
	"${TrimmedFile1}" \
	"${TrimmedFile2}" \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.sam

echo "Done."



echo -e "\n\n"
echo "###############################################"
echo "#Analysis step 3:  Marking PCR duplicate reads#"
echo "###############################################"
echo "Beginning analysis step 3 (marking and removing PCR duplicates) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

samtools view -bS \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.sam \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.bam

samtools sort \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.bam \
	-o "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted.bam

picard MarkDuplicates \
	INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted.bam \
	OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.bam \
	ASSUME_SORTED=true \
	METRICS_FILE="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.txt \
	VALIDATION_STRINGENCY=SILENT \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

samtools index \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.bam

samtools view -b \
	-F 0x400 \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.bam \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.bam

samtools index \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.bam

echo "Done."



echo -e "\n\n"
echo "####################################################"
echo "#Analysis step 4:  Remove low mapping quality reads#"
echo "####################################################"
echo "Beginning analysis step 4 (remove low mapping quality reads) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

samtools view \
	-bq "${FilterMappingQualityCutoff}" \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.bam \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam

samtools index \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam

echo "Done."



echo -e "\n\n"
echo "######################"
echo "#Analysis step 5:  QC#"
echo "######################"
echo "Beginning analysis step 5 (QC) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

samtools flagstat \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Alignment.txt

bedtools intersect \
	-v \
	-bed \
	-abam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-b "${BEDorGFFfile}" | \
	wc -l \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_ReadsOffTarget.txt

#Output depth of coverage for all regions in the BAM file
#Sequential positions at the same read depth are merged into a single region
bedtools genomecov \
	-ibam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-bga -split \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_CoverageTotal.bedgraph

#Output the mean depth of coverage for each region in the BED/GFF file
bedtools coverage \
	-a "${BEDorGFFfile}" \
	-b "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-mean \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_MeanCoverageBED.bedgraph

#Get percent genome covered
zero=$(bedtools genomecov \
	-ibam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-g "${Ref_FASTA}" \
	-bga | \
	awk '$4==0 {bpCountZero+=($3-$2)} {print bpCountZero}' | \
	tail -1)
nonzero=$(bedtools genomecov \
	-ibam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-g "${Ref_FASTA}" \
	-bga | \
	awk '$4>0 {bpCountNonZero+=($3-$2)} {print bpCountNonZero}' | \
	tail -1)
percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
echo -e "Number of bases at 0 read-depth:\t""${zero}" \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt
echo -e "Number of bases at >0 read-depth:\t""${nonzero}" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt
echo -e "Percent reference genome covered:\t""${percent}" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt

echo "Done."



echo -e "\n\n"
echo "#######################################################"
echo "#Analysis step 6:  Downsampling / random read sampling#"
echo "#######################################################"
echo "Beginning analysis step 6 (downsampling / random read sampling) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

picard DownsampleSam \
	INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam \
	RANDOM_SEED=50 \
	PROBABILITY=0.75 \
	VALIDATION_STRINGENCY=SILENT

samtools index \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam

picard DownsampleSam \
	INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam \
	RANDOM_SEED=50 \
	PROBABILITY=0.5 \
	VALIDATION_STRINGENCY=SILENT

samtools index \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam

picard DownsampleSam \
	INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam \
	RANDOM_SEED=50 \
	PROBABILITY=0.25 \
	VALIDATION_STRINGENCY=SILENT

samtools index \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam

echo "Done."



echo -e "\n\n"
echo "###################################"
echo "#Analysis step 7:  Variant calling#"
echo "###################################"
echo "Beginning analysis step 7 (variant calling) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

if [ $CallingStringency == "relaxed" ]
then
	bcftools mpileup \
	--redo-BAQ \
	--min-BQ 30 \
	--per-sample-mF \
	--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
	-f "${Ref_FASTA}" \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam | \
	bcftools call \
		--consensus-caller \
		--variants-only \
		--pval-threshold 1.0 \
		-Ob \
		> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf
elif [ $CallingStringency == "normal" ]
then
	bcftools mpileup \
		--redo-BAQ \
		--min-BQ 30 \
		--per-sample-mF \
		--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
		-f "${Ref_FASTA}" \
		"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
		"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam \
		"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam \
		"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam | \
		bcftools call \
			--multiallelic-caller \
			--variants-only \
			-Ob \
			> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf
fi

#Re-header the VCF
echo -e ""${PatientID}"" \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.txt
echo -e ""${PatientID}"_75pc" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.txt
echo -e ""${PatientID}"_50pc" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.txt
echo -e ""${PatientID}"_25pc" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.txt
bcftools reheader \
	--samples "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.txt \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf.tmp

cp -f "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf.tmp \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf

# 1st pipe, bcftools norm -Ou -m -any
#	Left-align and normalize indels; split multiallelic sites into multiple rows
#	-Ou, output in uncompressed format
#	-m, for multialleles
#	-any, SNPs and indels should be merged into a single record
#
# 2nd pipe, bcftools norm -Ou -f
#	As above but will check if reference alleles match the reference
#	-Ou, output in uncompressed format
#	-f, specify reference FASTA
bcftools norm \
	-Ou \
	-m \
	-any \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf \
	| bcftools norm \
		-Ov \
		-f "${Ref_FASTA}" \
		> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.vcf

#Filter the variants further
#Options:
#	-Q INT	minimum RMS mapping quality for SNPs [10]
#	-d INT	minimum read depth [2]
#	-D INT	maximum read depth [10000000]
#	-a INT	minimum number of alternate bases [2]
#	-w INT	SNP within INT bp around a gap to be filtered [3]
#	-W INT	window size for filtering adjacent gaps [10]
#	-1 FLOAT	min P-value for strand bias (given PV4) [0.0001]
#	-2 FLOAT	min P-value for baseQ bias [1e-100]
#	-3 FLOAT	min P-value for mapQ bias [0]
#	-4 FLOAT	min P-value for end distance bias [0.0001]
#	-e FLOAT	min P-value for HWE (plus F<0) [0.0001]
#	-p		print filtered variants
echo "Applying further filtering to called variants..."
echo "Variants filtered out:"
bcftools view \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.vcf \
	| vcfutils.pl varFilter \
		-d "${FilterReadDepthCutoff}" \
		-a 1 \
		-w 1 \
		-W 3 \
		-1 0.05 \
		-2 0.05 \
		-3 0.05 \
		-4 0.05 \
		-e 0.05 \
		-p \
		> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_FiltExtra.vcf

#Sort the VCF
bcftools sort -Ov \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_FiltExtra.vcf \
	> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Final.vcf

echo "Done."



echo -e "\n\n"
echo "##############################"
echo "#Analysis step 8:  Annotation#"
echo "##############################"
echo "Beginning analysis step 8 (annotation) on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

vep \
	--input_file "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Final.vcf \
	--format vcf \
	--species homo_sapiens \
	--assembly GRCh38 \
	--cache \
	--refseq \
	--offline \
	--fasta "${Ref_FASTA}" \
	--check_ref \
	--bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	--output_file "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.txt \
	--tab \
	--stats_file "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.html \
	--variant_class \
	--sift b \
	--polyphen b \
	--nearest symbol \
	--gene_phenotype \
	--regulatory \
	--numbers \
	--domains \
	--vcf_info_field VEP \
	--hgvs \
	--hgvsg \
	--symbol \
	--tsl \
	--canonical \
	--af \
	--af_1kg \
	--af_esp \
	--af_gnomad


echo -e "\n\n\n"
echo "############################################################"
echo "#Post-analysis step:  Moving files from temporary directory#"
echo "############################################################"
echo "Beginning post-analysis tidy-up on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

#Copy analysis files
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam.bai \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.txt \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PCRDuplicates.txt
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Alignment.txt \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_ReadsOffTarget.txt \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_CoverageTotal.bedgraph \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_MeanCoverageBED.bedgraph \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Final.vcf \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.txt \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.html \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${TrimmingReport1}" \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${TrimmingReport2}" \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${TrimmingHtml1}" \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"
mv "${TrimmingHtml2}" \
	"${Results_root}"/"${RunNumber}"/"${PatientID}"

rm -R "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/

echo "Done."

echo "Analysis script finished on `date`" \
	>> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PipelineLog.txt

exit 0