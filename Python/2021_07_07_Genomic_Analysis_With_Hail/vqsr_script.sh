#!/usr/bin/env bash

gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
    -V cohort.vcf.gz \
    --filter-expression "ExcessHet > 54.69" \
    --filter-name ExcessHet \
    -O cohort_excesshet.vcf.gz 

gatk MakeSitesOnlyVcf \
    -I cohort_excesshet.vcf.gz \
    -O cohort_sitesonly.vcf.gz

gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
    -V cohort_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    --mode INDEL \
    --max-gaussians 4 \
    -resource:mills,known=false,training=true,truth=true,prior=12 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -O cohort_indels.recal \
    --tranches-file cohort_indels.tranches

gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
    -V cohort_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -mode SNP \
    --max-gaussians 6 \
    -resource:hapmap,known=false,training=true,truth=true,prior=15 hapmap_3.3.hg38.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12 1000G_omni2.5.hg38.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -O cohort_snps.recal \
    --tranches-file cohort_snps.tranches

gatk --java-options "-Xmx5g -Xms5g" \
    ApplyVQSR \
    -V test_output.vcf.gz \
    --recal-file cohort_indels.recal \
    --tranches-file cohort_indels.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL \
    -O indel.recalibrated.vcf.gz

gatk --java-options "-Xmx5g -Xms5g" \
    ApplyVQSR \
    -V indel.recalibrated.vcf.gz \
    --recal-file cohort_snps.recal \
    --tranches-file cohort_snps.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode SNP \
    -O snp.recalibrated.vcf.gz