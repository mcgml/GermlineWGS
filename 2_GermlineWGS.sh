#!/bin/bash
#PBS -l walltime=999:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, Harriet Jackson All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="1.3.0"

# Script 2 runs in panel folder, requires final Bams, gVCFs and a PED file
# Variant filtering assumes non-related samples. If familiy structures are known they MUST be provided in the PED file

#load run & pipeline variables
. variables

addMetaDataToVCF(){
    output=$(echo "$1" | sed 's/\.vcf/_meta\.vcf/g')
    grep '^##' "$1" > "$output"
    for sample in $(/share/apps/bcftools-distros/bcftools-1.4/bcftools query -l "$1"); do
        cat "$sample"/"$seqId"_"$sample"_meta.txt >> "$output"
    done
    grep -v '^##' "$1" >> "$output"
}

annotateVCF(){
    #annotate VCF
    perl /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl \
    --verbose \
    --no_progress \
    --everything \
    --fork 12 \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file "$1" \
    --format vcf \
    --output_file "$2" \
    --force_overwrite \
    --no_stats \
    --cache \
    --dir /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --fasta /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --no_intergenic \
    --offline \
    --cache_version 86 \
    --allele_number \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --refseq
    
    #check VEP has produced annotated VCF
    if [ ! -e "$2" ]; then
        cp "$1" "$2"
    fi

    #index annotated VCF
    /share/apps/igvtools-distros/igvtools_2.3.75/igvtools index "$2"
}

### Joint variant calling and filtering ###

#Joint genotyping
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V GVCFs.list \
-o "$seqId"_variants.vcf \
-ped "$seqId"_pedigree.ped \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-nt 12

#Build the SNP recalibration model
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_variants.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /state/partition1/db/human/gatk/2.8/b37/hapmap_3.3.b37.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 /state/partition1/db/human/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile "$seqId"_SNP.recal \
-tranchesFile "$seqId"_SNP.tranches \
-rscriptFile "$seqId"_SNP_plots.R \
-ped "$seqId"_pedigree.ped \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-nt 12

#Apply the desired level of recalibration to the SNPs in the call set
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_variants.vcf \
-mode SNP \
--ts_filter_level 99.5 \
-recalFile "$seqId"_SNP.recal \
-tranchesFile "$seqId"_SNP.tranches \
-o "$seqId"_recalibrated_snps_raw_indels.vcf \
-ped "$seqId"_pedigree.ped \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-nt 12

#Build the Indel recalibration model
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_recalibrated_snps_raw_indels.vcf \
-resource:mills,known=false,training=true,truth=true,prior=12.0 /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--maxGaussians 4 \
-recalFile "$seqId"_INDEL.recal \
-tranchesFile "$seqId"_INDEL.tranches \
-rscriptFile "$seqId"_INDEL_plots.R \
-ped "$seqId"_pedigree.ped \
-nt 12

#Apply the desired level of recalibration to the Indels in the call set
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_recalibrated_snps_raw_indels.vcf \
-mode INDEL \
--ts_filter_level 99.5 \
-recalFile "$seqId"_INDEL.recal \
-tranchesFile "$seqId"_INDEL.tranches \
-o "$seqId"_recalibrated_variants.vcf \
-ped "$seqId"_pedigree.ped \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-nt 12

#Apply only family priors to a callset
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T CalculateGenotypePosteriors \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_recalibrated_variants.vcf \
--skipPopulationPriors \
-ped "$seqId"_pedigree.ped \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-o "$seqId"_recalibrated_variants_gcp.vcf

#phase by transmission
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T PhaseByTransmission \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_recalibrated_variants_gcp.vcf \
-ped "$seqId"_pedigree.ped \
-o "$seqId"_recalibrated_variants_gcp_phased.vcf \
--DeNovoPrior 0.000001 \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-mvf "$seqId"_MendelianViolations.txt

#filter genotypes
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_recalibrated_variants_gcp_phased.vcf \
-ped "$seqId"_pedigree.ped \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
--setFilteredGtToNocall \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-o "$seqId"_recalibrated_variants_gcp_phased_gtfiltered.vcf

#clean up
#TODO