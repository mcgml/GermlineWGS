#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, Harriet Jackson All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.0.0"

# Script 1 runs in sample folder, requires fastq files split by lane

#load sample & pipeline variables
. *.variables

countQCFlagFails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

### Preprocessing ###

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls *.fastq.gz | cut -d_ -f1-2 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f2)
    read1Fastq=$(ls "$fastqPair"_1.fastq.gz)
    read2Fastq=$(ls "$fastqPair"_2.fastq.gz)

    #trim adapters
    #TODO reaplace trimming with SparkTrim tool?
    /share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
    -m 35 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    "$read1Fastq" \
    "$read2Fastq"

    #FASTQ screen for inter-species contamination
    /share/apps/fastqscreen-distros/fastq_screen_v0.10.0/fastq_screen \
    --aligner bwa \
    --threads 12 \
    "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #convert fastq to ubam
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar FastqToSam \
    F1="$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    F2="$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    O="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
    READ_GROUP_NAME="$seqId"_"$laneId"_"$sampleId" \
    SAMPLE_NAME="$sampleId" \
    LIBRARY_NAME="$worklistId"_"$sampleId" \
    PLATFORM_UNIT="$seqId"_"$laneId" \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="WTCHG" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/state/partition1/tmpdir

    #fastqc
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq *_fastqc.zip *_fastqc.html

done

#merge lane bams
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MergeSamFiles \
$(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
SORT_ORDER=queryname \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
USE_THREADING=true \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir \
O="$seqId"_"$sampleId"_unaligned.bam

#uBam2fq, map & MergeBamAlignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
INTERLEAVE=true \
NON_PF=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
COMPRESSION_LEVEL=0 \
TMP_DIR=/state/partition1/tmpdir | \
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-t 12 \
-p \
/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37_decoy_phix.fasta \
/dev/stdin | \
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MergeBamAlignment \
EXPECTED_ORIENTATIONS=FR \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \
OUTPUT="$seqId"_"$sampleId"_aligned.bam \
REFERENCE_SEQUENCE=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37_decoy_phix.fasta \
PAIRED_RUN=true \
SORT_ORDER="coordinate" \
IS_BISULFITE_SEQUENCE=false \
ALIGNED_READS_ONLY=false \
CLIP_ADAPTERS=false \
MAX_RECORDS_IN_RAM=2000000 \
ADD_MATE_CIGAR=true \
MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
UNMAP_CONTAMINANT_READS=false \
CLIP_OVERLAPPING_READS=true \
ALIGNER_PROPER_PAIR_FLAGS=false \
INCLUDE_SECONDARY_ALIGNMENTS=true \
CREATE_INDEX=true \
TMP_DIR=/state/partition1/tmpdir

if [ "$pcr" = true ]; then
    #Mark duplicate reads
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MarkDuplicates \
    INPUT="$seqId"_"$sampleId"_aligned.bam \
    OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
    METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
    CREATE_INDEX=true \
    MAX_RECORDS_IN_RAM=2000000 \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=/state/partition1/tmpdir
else 
    mv "$seqId"_"$sampleId"_aligned.bam "$seqId"_"$sampleId"_rmdup.bam
    mv "$seqId"_"$sampleId"_aligned.bai "$seqId"_"$sampleId"_rmdup.bai
fi

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx24g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realign.intervals \
-nt 12 \
-XL NC_007605 -XL hs37d5 -XL phix \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-LOD 0.4 \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_realign.intervals \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realigned.bam \
-dt NONE

#Analyse patterns of covariation in the sequence dataset
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx6g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_realigned.bam \
-o "$seqId"_"$sampleId"_recal_data.table \
-L 22 \
-nct 12 \
-dt NONE

#Do a second pass to analyze covariation remaining after recalibration
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx6g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-BQSR "$seqId"_"$sampleId"_recal_data.table \
-I "$seqId"_"$sampleId"_realigned.bam \
-o "$seqId"_"$sampleId"_post_recal_data.table \
-L 22 \
-nct 12 \
-dt NONE

#Generate BQSR plots
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-before "$seqId"_"$sampleId"_recal_data.table \
-after "$seqId"_"$sampleId"_post_recal_data.table \
-plots "$seqId"_"$sampleId"_recalibration_plots.pdf \
-csv "$seqId"_"$sampleId"_recalibration.csv \
-dt NONE

#Apply the recalibration to your sequence data
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-I "$seqId"_"$sampleId"_realigned.bam \
-BQSR "$seqId"_"$sampleId"_recal_data.table \
-o "$seqId"_"$sampleId".bam \
-nct 6 \
-dt NONE

### Variant calling ###

#SNPs and Indels GVCF with Haplotypecaller
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-I "$seqId"_"$sampleId".bam \
-o "$seqId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
--emitRefConfidence GVCF \
-nct 12 \
-XL NC_007605 -XL hs37d5 -XL phix \
-dt NONE

### QC ###

#Alignment metrics: library sequence similarity
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectAlignmentSummaryMetrics \
R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_AlignmentSummaryMetrics.txt \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir

#Calculate insert size: fragmentation performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectInsertSizeMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_InsertMetrics.txt \
H="$seqId"_"$sampleId"_InsertMetrics.pdf \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir

#Generate per-base coverage: variant detection sensitivity
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx12g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37_decoy_phix.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
--minBaseQuality 10 \
-ct 20 \
--omitIntervalStatistics \
--omitLocusTable \
-rf MappingQualityUnavailable \
-nt 12 \
-XL NC_007605 -XL hs37d5 -XL phix \
-dt NONE

#tabix index the per-base coverage file
awk -F'[\t|:]' '{if(NR>1) print $1"\t"$2"\t"$3}' "$seqId"_"$sampleId"_DepthOfCoverage | \
/share/apps/htslib-distros/htslib-1.4/bgzip > "$seqId"_"$sampleId"_DepthOfCoverage.gz
/share/apps/htslib-distros/htslib-1.4/tabix -b2 -e2 -s1 "$seqId"_"$sampleId"_DepthOfCoverage.gz

#calculate CDS coverage with padding
zcat /state/partition1/db/human/refseq/ref_GRCh37.p13_top_level_canonical_b37_sorted.gff3.gz | \
grep "NP_[0-9]*\.[0-9]*" | \
awk -F'[\t|;|=]' -v p=5 '$3 == "CDS" { gene="null"; for (i=9;i<NF;i++) if ($i=="gene"){gene=$(i+1); break}; genes[gene] = genes[gene]$1"\t"($4-1)-p"\t"$5+p"\t"gene";" } END { for (gene in genes) print genes[gene] }' | \
while read line; do
    echo "$line" | \
    tr ';' '\n'| \
    sort -k1,1V -k2,2n -k3,3n | \
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge -c 4 -o distinct;
done | \
sort -k1,1V -k2,2n -k3,3n > "$panel"_ClinicalCoverageTargets.bed

#Make PASS BED
/share/apps/htslib-distros/htslib-1.4/tabix -R "$panel"_ClinicalCoverageTargets.bed \
"$seqId"_"$sampleId"_DepthOfCoverage.gz | \
awk -v '$3 >= 20 { print $1"\t"$2-1"\t"$2 }' | \
sort -k1,1V -k2,2n -k3,3n | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge > "$seqId"_"$sampleId"_PASS.bed

#Calculate overlap between PASS BED and ClinicalCoverageTargets
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools coverage \
-a "$panel"_ClinicalCoverageTargets.bed \
-b "$seqId"_"$sampleId"_PASS.bed | \
awk '{pass[$4]+=$6; len[$4]+=$7} END { for(i in pass) printf "%s\t %.2f%\n", i, (pass[i]/len[i]) * 100 }' | \
sort -k1,1 > "$seqId"_"$sampleId"_ClinicalCoverageGeneCoverage.txt

#Calculate dna contamination: sample-to-sample contamination
#TODO check args
/share/apps/verifyBamID-distros/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID \
--vcf /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
--bam "$seqId"_"$sampleId".bam \
--out "$seqId"_"$sampleId"_Contamination \
--verbose \
--ignoreRG \
--chip-none \
--minMapQ 20 \
--maxDepth 50 \
--precise

#create final file lists
find $PWD -name "$seqId"_"$sampleId".g.vcf >> ../GVCFs.list
find $PWD -name "$seqId"_"$sampleId".bam >> ../BAMs.list

#check if all VCFs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../GVCFs.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo -e "seqId=$seqId\npanel=$panel" > ../variables
    cp 2_GermlineWGS.sh .. && cd .. && qsub 2_GermlineWGS.sh
fi

#clean up
#TODO