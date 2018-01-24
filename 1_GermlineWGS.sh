#!/bin/bash
#PBS -l walltime=999:00:00
#PBS -l ncpus=12
. /share/apps/root-distros/root_v5.34.36.Linux-slc6-x86_64-gcc4.4/bin/thisroot.sh
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, Harriet Jackson All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.3.0"

# Script 1 runs in sample folder, requires fastq files split by lane
# Designed for a single germline WGS analysis from 1 lane of data ~ 120Gb

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
setOpticalDuplicatePixelDistance() {
    #determine the pixel distance depending on the flowcell
    if [ "$1" = "true" ]; then
        echo -n "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
    else
        echo -n "OPTICAL_DUPLICATE_PIXEL_DISTANCE=100"
    fi
}
setPcrIndelModel(){
    #determine the indel model depending on the library type
    if [ "$1" = "true" ]; then
        echo -n "--pcr_indel_model CONSERVATIVE"
    else
        echo -n "--pcr_indel_model NONE"
    fi
}

### Preprocessing ###

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

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
    SEQUENCING_CENTER="$sequencingCentre" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/state/partition1/tmpdir

    #uBam2fq, map & MergeBamAlignment
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar SamToFastq \
    I="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
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
    /state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
    /dev/stdin | \
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MergeBamAlignment \
    EXPECTED_ORIENTATIONS=FR \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
    OUTPUT="$seqId"_"$sampleId"_"$laneId"_aligned.bam \
    REFERENCE_SEQUENCE=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
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

    #fastqc
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi

    #FASTQ screen for inter-species contamination
    /share/apps/fastqscreen-distros/fastq_screen_v0.10.0/fastq_screen \
    --aligner bwa \
    --threads 12 \
    "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq *_fastqc.zip *_fastqc.html

done

#Mark duplicate reads & merge aBAMs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MarkDuplicates \
$(ls "$seqId"_"$sampleId"_*_aligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
$(setOpticalDuplicatePixelDistance "$patterned") \
TMP_DIR=/state/partition1/tmpdir

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realign.intervals \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-nt 12

#Realign around indels
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_realign.intervals \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realigned.bam

#Analyse patterns of covariation in the sequence dataset
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx6g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_realigned.bam \
-o "$seqId"_"$sampleId"_recal_data.table \
-L 20 \
-nct 12

#Do a second pass to analyze covariation remaining after recalibration
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx6g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-BQSR "$seqId"_"$sampleId"_recal_data.table \
-I "$seqId"_"$sampleId"_realigned.bam \
-o "$seqId"_"$sampleId"_post_recal_data.table \
-L 20 \
-nct 12

#Generate BQSR plots
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-before "$seqId"_"$sampleId"_recal_data.table \
-after "$seqId"_"$sampleId"_post_recal_data.table \
-plots "$seqId"_"$sampleId"_recalibration_plots.pdf \
-csv "$seqId"_"$sampleId"_recalibration.csv

#Apply the recalibration to your sequence data
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_realigned.bam \
-BQSR "$seqId"_"$sampleId"_recal_data.table \
-o "$seqId"_"$sampleId".bam \
-nct 6

### Variant calling ###

#SNPs and Indels GVCF with Haplotypecaller
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx32g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId".bam \
-o "$seqId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
--emitRefConfidence GVCF \
$(setPcrIndelModel "$pcr") \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-nct 12

#Structural variant calling with Manta
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools subtract \
-a <(awk -F"\t" '{print $1"\t"0"\t"$2}' /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai) \
-b /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed | \
/share/apps/htslib-distros/htslib-1.4.1/bgzip -c > manta.bed.gz
/share/apps/htslib-distros/htslib-1.4.1/tabix -p bed manta.bed.gz

/share/apps/manta-distros/manta-1.2.1.centos6_x86_64/bin/configManta.py \
--bam "$seqId"_"$sampleId".bam \
--referenceFasta /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--callRegions manta.bed.gz \
--runDir manta
manta/runWorkflow.py \
--quiet \
-m local \
-j 12

#soft link VCF
ln -s manta/results/variants/diploidSV.vcf.gz "$seqId"_"$sampleId"_sv_filtered.vcf.gz
ln -s manta/results/variants/diploidSV.vcf.gz.tbi "$seqId"_"$sampleId"_sv_filtered.vcf.gz.tbi

#RD-CNV calling with CNVNator

#extract reads
/share/apps/CNVnator-distros/CNVnator_v0.3.3/src/cnvnator \
-chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
-root "$seqId"_"$sampleId".root \
-tree "$seqId"_"$sampleId".bam \
-unique

#create histogram
/share/apps/CNVnator-distros/CNVnator_v0.3.3/src/cnvnator \
-chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
-root "$seqId"_"$sampleId".root \
-his 100 \
-d /share/apps/CNVnator-distros/gatk_v0_hg38_split

#stats
/share/apps/CNVnator-distros/CNVnator_v0.3.3/src/cnvnator \
-chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
-root "$seqId"_"$sampleId".root \
-stat 100

#partition
/share/apps/CNVnator-distros/CNVnator_v0.3.3/src/cnvnator \
-chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
-root "$seqId"_"$sampleId".root \
-partition 100

#calling
echo -e "CNV_type\tcoordinates\tCNV_size\tnormalized_RD\te-val1\te-val2\te-val3\te-val4\tq0" > "$seqId"_"$sampleId"_cnv.txt
/share/apps/CNVnator-distros/CNVnator_v0.3.3/src/cnvnator \
-chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
-root "$seqId"_"$sampleId".root \
-call 100 >> "$seqId"_"$sampleId"_cnv.txt

### QC ###

#Alignment metrics: library sequence similarity
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectAlignmentSummaryMetrics \
R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
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

#Coverage analysis
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectWgsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_WgsMetrics.txt \
MINIMUM_MAPPING_QUALITY=20 \
MINIMUM_BASE_QUALITY=10 \
R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta

#find callable regions
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T CallableLoci \
--maxDepth 65 \
--minDepth 10 \
--minBaseQuality 10 \
--minMappingQuality 20 \
--format BED \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId".bam \
-o "$seqId"_"$sampleId"_callable_status.bed \
--summary "$seqId"_"$sampleId"_callable_status.txt \
-XL /data/diagnostics/pipelines/GermlineWGS/GermlineWGS-"$version"/not_callable.bed \
-rf MappingQualityUnavailable

#make callable bed
awk -F"\t" '$4=="CALLABLE"' "$seqId"_"$sampleId"_callable_status.bed | \
/share/apps/htslib-distros/htslib-1.4.1/bgzip -c > "$seqId"_"$sampleId"_callable.bed.gz
/share/apps/htslib-distros/htslib-1.4.1/tabix -p bed "$seqId"_"$sampleId"_callable.bed.gz

### Clean up ###

#create final file lists
find $PWD -name "$seqId"_"$sampleId".g.vcf >> ../GVCFs.list
find $PWD -name "$seqId"_"$sampleId".bam >> ../BAMs.list

#check if all VCFs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../GVCFs.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo "seqId=$seqId" > ../variables
    cp 2_GermlineWGS.sh .. && cd .. && qsub 2_GermlineWGS.sh
fi

#delete unnecessary files
#TODO