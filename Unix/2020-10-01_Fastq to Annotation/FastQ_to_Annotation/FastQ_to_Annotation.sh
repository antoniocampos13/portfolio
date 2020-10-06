#!/usr/bin/env bash
# GPL License and Copyright Notice #######################################################
# Copyright (C)  2019  Kevin Blighe: Original implementation
# Copyright (C)  2020  Antonio Victor Campos Coelho: Modified version.
#
# Permission is granted to copy, distribute and/or modify this document
# under the terms of the GNU Free Documentation License, Version 1.3
# or any later version published by the Free Software Foundation;
# with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
# A copy of the license is included in the section entitled "GNU
# Free Documentation License".
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################################

#########################################################
#Pre-liminary step:  Safeguarding the script's operation#
#########################################################
echo "Checking command line parameters..."

Read1="${1}_1.fastq.gz"
Read2="${1}_2.fastq.gz"
Ref_FASTA="${2}"
BEDorGFFfile="${3}"
TrimmingQualityReadEnds="${4}"
TrimmingReadLengthMin="${5}"
TrimmingAdaptor="${6}"
FilterReadDepthCutoff="${7}"
FilterMappingQualityCutOff="${8}"
CallingStringency="${9}"
User="${10}"


#Check that 10 parameters have been passed
if [ $# -ne 10 ]
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

echo "Running script ..."

echo "Log initialisation on `date`..."
echo "Log file is "${1}"_Master_Log.txt"

#Begin log file
echo "Beginning script on `date`, run by "${User}" with the following parameters:" \
	> "${1}"_Pipeline_Log.txt
echo -e "\t1\t"${Read1}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t2\t"${Read2}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t3\t"${Ref_FASTA}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t4\t"${RunNumber}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t5\t"${PatientID}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t6\t"${BEDorGFFfile}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t7\t"${TrimmingQualityReadEnds}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t8\t"${TrimmingReadLengthMin}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t9\t"${TrimmingAdaptor}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t10\t"${FilterReadDepthCutoff}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t11\t"${FilterMappingQualityCutOff}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t12\t"${CallingStringency}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t13\t"${Results_root}"" \
	>> "${1}"_Pipeline_Log.txt
echo -e "\t14\t"${User}"" \
	>> "${1}"_Pipeline_Log.txt

exec &> "${1}"_Master_Log.txt

echo -e "\n\n"
echo "#####################################################"
echo "#Step 1:  adaptor and read quality trimming         #"
echo "#####################################################"
echo "Beginning Step 1 (adaptor and read quality trimming) on `date`" \
	>> "${1}"_Pipeline_Log.txt

trim_galore \
	--qual "${TrimmingQualityReadEnds}" \
	--length "${TrimmingReadLengthMin}" \
	--"${TrimmingAdaptor}" \
	--paired "${Read1}" "${Read2}" \
	--fastqc

TrimmedFile1="${1}"_1_val_1.fq.gz
TrimmedFile2="${1}"_2_val_2.fq.gz

echo "Done."
echo -e "\n\n"
echo "#############################"
echo "#Step 2:  Alignment         #"
echo "#############################"
echo "Beginning Step 2 (alignment) on `date`" \
	>> "${1}"_Pipeline_Log.txt

#Producing aligned SAM files with BWA mem
bwa mem \
	"${Ref_FASTA}" \
	"${TrimmedFile1}" \
	"${TrimmedFile2}" \
	> "${1}"_Aligned.sam

echo "Done."

echo -e "\n\n"
echo "###############################################"
echo "#Step 3:  Marking PCR duplicate reads         #"
echo "###############################################"
echo "Beginning Step 3 (marking and removing PCR duplicates) on `date`" \
	>> "${1}"_Pipeline_Log.txt

samtools view -bS \
	"${1}"_Aligned.sam \
	> "${1}"_Aligned.bam

samtools sort \
	"${1}"_Aligned.bam \
	-o "${1}"_Aligned_Sorted.bam

picard MarkDuplicates \
	INPUT="${1}"_Aligned_Sorted.bam \
	OUTPUT="${1}"_Aligned_Sorted_PCRDupes.bam \
	ASSUME_SORTED=true \
	METRICS_FILE="${1}"_Aligned_Sorted_PCRDupes.txt \
	VALIDATION_STRINGENCY=SILENT \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

samtools index \
	"${1}"_Aligned_Sorted_PCRDupes.bam

samtools view -b \
	-F 0x400 \
	"${1}"_Aligned_Sorted_PCRDupes.bam \
	> "${1}"_Aligned_Sorted_PCRDuped.bam

samtools index \
	"${1}"_Aligned_Sorted_PCRDuped.bam

echo "Done."

echo -e "\n\n"
echo "####################################################"
echo "#Step 4:  Remove low mapping quality reads         #"
echo "####################################################"
echo "Beginning Step 4 (remove low mapping quality reads) on `date`" \
	>> "${1}"_Pipeline_Log.txt

samtools view \
	-bq "${FilterMappingQualityCutoff}" \
	"${1}"_Aligned_Sorted_PCRDuped.bam \
	> "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam

samtools index \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam

echo "Done."
echo -e "\n\n"
echo "######################"
echo "#Step 5:  QC         #"
echo "######################"
echo "Beginning Step 5 (QC) on `date`" \
	>> "${1}"_Pipeline_Log.txt

samtools flagstat \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	> "${1}"_Alignment.txt

bedtools intersect \
	-v \
	-bed \
	-abam "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-b "${BEDorGFFfile}" | \
	wc -l \
	> "${1}"_ReadsOffTarget.txt

#Output depth of coverage for all regions in the BAM file
#Sequential positions at the same read depth are merged into a single region
bedtools genomecov \
	-ibam "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-bga -split \
	> "${1}"_CoverageTotal.bedgraph

#Output the mean depth of coverage for each region in the BED/GFF file
bedtools coverage \
	-a "${BEDorGFFfile}" \
	-b "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-mean \
	> "${1}"_MeanCoverageBED.bedgraph

#Get percent genome covered
zero=$(bedtools genomecov \
	-ibam "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-g "${Ref_FASTA}" \
	-bga | \
	awk '$4==0 {bpCountZero+=($3-$2)} {print bpCountZero}' | \
	tail -1)
nonzero=$(bedtools genomecov \
	-ibam "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	-g "${Ref_FASTA}" \
	-bga | \
	awk '$4>0 {bpCountNonZero+=($3-$2)} {print bpCountNonZero}' | \
	tail -1)
percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
echo -e "Number of bases at 0 read-depth:\t""${zero}" \
	> "${1}"_PercentGenomeCovered.txt
echo -e "Number of bases at >0 read-depth:\t""${nonzero}" \
	>> "${1}"_PercentGenomeCovered.txt
echo -e "Percent reference genome covered:\t""${percent}" \
	>> "${1}"_PercentGenomeCovered.txt

echo "Done."

echo -e "\n\n"
echo "#######################################################"
echo "#Step 6:  Downsampling / random read sampling         #"
echo "#######################################################"
echo "Beginning Step 6 (downsampling / random read sampling) on `date`" \
	>> "${1}"_Pipeline_Log.txt

picard DownsampleSam \
	INPUT="${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	OUTPUT="${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam \
	RANDOM_SEED=50 \
	PROBABILITY=0.75 \
	VALIDATION_STRINGENCY=SILENT

samtools index \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam

picard DownsampleSam \
	INPUT="${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	OUTPUT="${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam \
	RANDOM_SEED=50 \
	PROBABILITY=0.5 \
	VALIDATION_STRINGENCY=SILENT

samtools index \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam

picard DownsampleSam \
	INPUT="${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	OUTPUT="${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam \
	RANDOM_SEED=50 \
	PROBABILITY=0.25 \
	VALIDATION_STRINGENCY=SILENT

samtools index \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam

echo "Done."

echo -e "\n\n"
echo "###################################"
echo "#Step 7:  Variant calling         #"
echo "###################################"
echo "Beginning Step 7 (variant calling) on `date`" \
	>> "${1}"_Pipeline_Log.txt

if [ $CallingStringency == "relaxed" ]
then
	bcftools mpileup \
	--redo-BAQ \
	--min-BQ 30 \
	--per-sample-mF \
	--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
	-f "${Ref_FASTA}" \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam | \
	bcftools call \
		--consensus-caller \
		--variants-only \
		--pval-threshold 1.0 \
		-Ob \
		> "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf
elif [ $CallingStringency == "normal" ]
then
	bcftools mpileup \
		--redo-BAQ \
		--min-BQ 30 \
		--per-sample-mF \
		--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
		-f "${Ref_FASTA}" \
		"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
		"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam \
		"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam \
		"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam | \
		bcftools call \
			--multiallelic-caller \
			--variants-only \
			-Ob \
			> "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf
fi

#Re-header the VCF
echo -e ""${1}"" \
	> samples.txt
echo -e ""${1}"_75pc" \
	>> samples.txt
echo -e ""${1}"_50pc" \
	>> samples.txt
echo -e ""${1}"_25pc" \
	>> samples.txt
bcftools reheader \
	--samples samples.txt \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf \
	> "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf.tmp

cp -f "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf.tmp \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf

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
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf \
	| bcftools norm \
		-Ov \
		-f "${Ref_FASTA}" \
		> "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.vcf

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
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.vcf \
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
		> "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_FiltExtra.vcf

#Sort the VCF
bcftools sort -Ov \
	"${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_FiltExtra.vcf \
	> "${1}"_Final.vcf

echo -e "\n\n"
echo "Done."

echo -e "\n\n"
echo "##############################"
echo "#Step 8:  Annotation         #"
echo "##############################"
echo "Beginning Step 8 (annotation) on `date`" \
	>> "${1}"_Pipeline_Log.txt

vep \
	--input_file "${1}"_Final.vcf \
	--format vcf \
	--species homo_sapiens \
	--assembly GRCh38 \
	--cache \
	--offline \
	--refseq \
	--fasta "${Ref_FASTA}" \
	--check_ref \
	--bam "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam \
	--output_file "${1}"_AnnotationVEP.txt \
	--tab \
	--stats_file "${1}"_AnnotationVEP.html \
	--everything

echo "Done."
echo "Script finished on `date`" \
	>> "${1}"_Pipeline_Log.txt

echo -e "\n\n\n"
echo "############################################################"
echo "#Tidy-up                                                   #"
echo "############################################################"
echo "Beginning tidy-up on `date`" \
	>> "${1}"_Pipeline_Log.txt

mkdir -p "${1}_results"
mkdir -p "${1}_results"/trimmed_files
mkdir -p "${1}_results"/alignment_files
mkdir -p "${1}_results"/variant_call_files

mv *_fastqc.html *_trimming_report.txt *.zip *.fq.gz "${1}_results"/trimmed_files

mv *.bam *.sam *.bcf *_Aligned_Sorted_PCRDupes.txt *.bai *.tmp *.bedgraph "${1}"_Alignment.txt "${1}"_PercentGenomeCovered.txt "${1}"_ReadsOffTarget.txt "${1}_results"/alignment_files

mv samples.txt "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ.vcf "${1}"_Aligned_Sorted_PCRDuped_FiltMAPQ_FiltExtra.vcf "${1}_results"/variant_call_files

mv "${1}"_Master_Log.txt "${1}"_Pipeline_Log.txt "${1}"_Final.vcf "${1}"_AnnotationVEP.txt "${1}"_AnnotationVEP.html "${1}_results"/

exit 0