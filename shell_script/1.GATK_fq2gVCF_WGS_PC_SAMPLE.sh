#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J 1.GATK_fq2gVCF_WGS_PC_SAMPLE
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92G
#SBATCH -o /home/sevemthday77/log/1.GATK_fq2gVCF_WGS_PC_SAMPLE.out.txt
#SBATCH -e /home/sevemthday77/log/1.GATK_fq2gVCF_WGS_PC_SAMPLE.err.txt
#SBATCH --mail-user=sevemthday@gmail.com
#SBATCH --mail-type=END,FAIL

# run info
sample=SAMPLE
date=DATE
base_path=/staging/biology/sevemthday77/PancreaticCancer/WGS_res

# path
fastq_path=$base_path/WGS_fastq/$date
bam_path=$base_path/WGS_bam_gvcf/$date
sv_call_path=$base_path/WGS_SV/$date

mkdir -p $bam_path

# alignment
fq1=$(ls $fastq_path/${sample}* |grep "R1")
fq2=$(ls $fastq_path/${sample}* |grep "R2")
readgroup="@RG\tID:${sample}\tSM:${sample}\tLB:WES\tPL:Illumina"
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
snp=/staging/biology/sevemthday77/GATK_bundle/dbsnp_146.hg38.vcf.gz
indelone=/staging/biology/sevemthday77/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
indeltwo=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz
bwa=/opt/ohpc/Taiwania3/pkg/biology/BWA/BWA_v0.7.17/bwa
samtools=/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools
Picard_path=/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.26.0
GATK=/opt/ohpc/Taiwania3/pkg/biology/GATK/gatk_v4.2.0.0/gatk
export PATH=/opt/ohpc/Taiwania3/pkg/biology/Python/Python_v2.7.18/bin:$PATH

cd $bam_path
mkdir $sample
cd $sample

$bwa mem -t 14 -M -R $readgroup $ref $fq1 $fq2 | $samtools view -hb -o ${sample}.bam -

rm -rf $fq1
rm -rf $fq2

java -Xmx77g -jar ${Picard_path}/picard.jar SortSam \
CREATE_INDEX=true \
INPUT=${sample}.bam \
OUTPUT=${sample}_sort.bam \
SORT_ORDER=coordinate

rm -rf ${sample}.bam

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./" MarkDuplicates \
-I ${sample}_sort.bam \
-O ${sample}_marked.bam \
-M $sample.metrics

rm -rf ${sample}_sort.bam
rm -rf ${sample}_sort.bai

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./ -XX:ParallelGCThreads=12" BaseRecalibrator \
-R $ref \
-I ${sample}_marked.bam \
--known-sites $snp \
--known-sites $indelone \
--known-sites $indeltwo \
-O ${sample}_recal.table

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./ -XX:ParallelGCThreads=12" ApplyBQSR \
-R $ref \
-I ${sample}_marked.bam \
-bqsr ${sample}_recal.table \
-O ${sample}_bqsr.bam

rm -rf ${sample}_marked.bam

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./ -XX:ParallelGCThreads=12" HaplotypeCaller \
-R $ref \
-I ${sample}_bqsr.bam \
--dbsnp $snp \
-O ${sample}_raw.g.vcf.gz \
-ERC GVCF


rm tmp_read_resource_*.config

# SV
bam=$bam_path/$sample/${sample}_bqsr.bam

lumpy_thread=14
gridss_thread=14
svaba_thread=14
cue_thread=14
fastqc_thread=14

mkdir -p $sv_call_path/manta/$sample
mkdir -p $sv_call_path/delly/$sample
mkdir -p $sv_call_path/lumpy/$sample
mkdir -p $sv_call_path/gridss/$sample
mkdir -p $sv_call_path/svaba/$sample
mkdir -p $sv_call_path/cue/$sample
mkdir -p $sv_call_path/survivor/$sample
mkdir -p $sv_call_path/annotsv/$sample

# manta
PATH_raw=$PATH
export PATH=/opt/ohpc/Taiwania3/pkg/biology/Python/Python_v2.7.18/bin:$PATH
ml old-module
ml biology/bcftools/1.13

ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
manta_folder_path=/opt/ohpc/Taiwania3/pkg/biology/Manta/Manta_v1.6.0

$manta_folder_path/bin/configManta.py \
--bam $bam \
--referenceFasta $ref \
--runDir $sv_call_path/manta/$sample

$sv_call_path/manta/$sample/runWorkflow.py

bcftools view -i 'FILTER="PASS"' $sv_call_path/manta/$sample/results/variants/diploidSV.vcf.gz > $sv_call_path/manta/$sample/${sample}_manta.pass.vcf

rm -rf $sv_call_path/manta/$sample/workspace

ml purge
export PATH=$PATH_raw

# delly
PATH_raw=$PATH
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
delly_exclude_region=/staging/biology/sevemthday77/software/delly/human.hg38.excl.tsv

ml old-module
ml biology/DELLY/1.1.5
ml biology/bcftools/1.13

delly call \
-g $ref \
-o $sv_call_path/delly/$sample/${sample}_delly.bcf \
-x $delly_exclude_region \
$bam

bcftools view $sv_call_path/delly/$sample/${sample}_delly.bcf -Ov -i'GT="alt" & FILTER="PASS"' -o $sv_call_path/delly/$sample/${sample}_delly.pass.vcf

ml purge
export PATH=$PATH_raw

# lumpy
PATH_raw=$PATH
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
lumpy_exclude_region=/staging/biology/sevemthday77/software/smoove/bed/exclude.cnvnator_100bp.GRCh38.20170403.bed
export PATH=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin:$PATH
export PATH=/staging/biology/sevemthday77/software/gsort/bin:$PATH
export PATH=/opt/ohpc/Taiwania3/pkg/biology/lumpy-sv/lumpy-sv_v0.3.1/bin:$PATH
export PATH=/staging/biology/sevemthday77/software/htslib/htslib-1.15.1/bin:$PATH
export PATH=/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin:$PATH
export PATH=/opt/ohpc/Taiwania3/pkg/biology/svtyper/svtyper_v0.7.1/bin:$PATH
export PATH=/staging/biology/sevemthday77/software/mosdepth/bin:$PATH
export PATH=/staging/biology/sevemthday77/software/duphold/bin:$PATH
export PATH=/staging/biology/sevemthday77/software/smoove/bin:$PATH

smoove call -x --genotype \
--name $sample \
--outdir $sv_call_path/lumpy/$sample \
-f $ref \
--processes $lumpy_thread \
--exclude $lumpy_exclude_region \
$bam

cp $sv_call_path/lumpy/$sample/$sample-smoove.genotyped.vcf.gz $sv_call_path/lumpy/$sample/${sample}_lumpy.vcf.gz
bgzip -d $sv_call_path/lumpy/$sample/${sample}_lumpy.vcf.gz

export PATH=$PATH_raw

rm $sv_call_path/lumpy/$sample/$sample.split.bam.orig.bam
rm $sv_call_path/lumpy/$sample/$sample.disc.bam.orig.bam
rm $sv_call_path/lumpy/$sample/$sample.split.bam
rm $sv_call_path/lumpy/$sample/$sample.split.bam.csi
rm $sv_call_path/lumpy/$sample/$sample.disc.bam
rm $sv_call_path/lumpy/$sample/$sample.disc.bam.csi

# gridss
PATH_raw=$PATH
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
gridss_exclude_region=/staging/biology/sevemthday77/software/gridss/bedfile/ENCFF356LFX.bed

ml old-module
ml biology/bcftools/1.13
source /opt/ohpc/Taiwania3/pkg/biology/GRIDSS/GRIDSS_v2.13.1/env.sh

current_path=$(pwd)

gridss \
-r $ref \
-j /opt/ohpc/Taiwania3/pkg/biology/GRIDSS/GRIDSS_v2.13.1/gridss-2.13.1-gridss-jar-with-dependencies.jar \
-o $sv_call_path/gridss/$sample/${sample}_gridss.vcf \
-b $gridss_exclude_region \
-t $gridss_thread \
$bam

bcftools view -i 'FILTER="PASS"' $sv_call_path/gridss/$sample/${sample}_gridss.vcf > $sv_call_path/gridss/$sample/${sample}_gridss.pass.vcf

ml R/4.3.3
simpleEventAnnotation=/staging/biology/sevemthday77/software/Rscript_/simple-event-annotation_T3_v2.R
Rscript $simpleEventAnnotation \
-i $sv_call_path/gridss/$sample/${sample}_gridss.pass.vcf \
-o $sv_call_path/gridss/$sample/${sample}_gridss.pass.SEA.vcf \
-g hg38

ml purge
export PATH=$PATH_raw

rm $sv_call_path/gridss/$sample/${sample}_gridss.vcf.assembly.bam

rm $current_path/${sample}_bqsr.bam.gridss.working/${sample}_bqsr.bam.sv.bam
rm $current_path/${sample}_bqsr.bam.gridss.working/${sample}_bqsr.bam.sv.bam.csi
mv $current_path/${sample}_bqsr.bam.gridss.working/ $sv_call_path/gridss/$sample/
rm -rf $current_path/${sample}_gridss.vcf.gridss.working
rm $current_path/${sample}_gridss.vcf.assembly.bam.gridss.working/${sample}_gridss.vcf.assembly.bam.sv.bam
rm $current_path/${sample}_gridss.vcf.assembly.bam.gridss.working/${sample}_gridss.vcf.assembly.bam.sv.bam.bai
mv $current_path/${sample}_gridss.vcf.assembly.bam.gridss.working/ $sv_call_path/gridss/$sample/
mv $current_path/gridss.full.*.log $sv_call_path/gridss/$sample/
mv $current_path/gridss.timing.*.log $sv_call_path/gridss/$sample/
rm $current_path/libsswjni.so
rm $current_path/${sample}_bqsr.bam.*.750.SHORT.auxindex



# svaba
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
svaba=/opt/ohpc/Taiwania3/pkg/biology/svaba/svaba_v1.1.0/bin/svaba
ml old-module
ml biology/svaba/1.1.0
ml biology/HTSLIB/1.18
ml biology/bcftools/1.13

svaba run \
-t $bam \
-p $svaba_thread \
-a $sv_call_path/svaba/$sample/$sample \
-L 6 \
-I \
-G $ref

bgzip -c $sv_call_path/svaba/$sample/$sample.svaba.sv.vcf > $sv_call_path/svaba/$sample/${sample}_svaba.vcf.gz
tabix -p vcf $sv_call_path/svaba/$sample/${sample}_svaba.vcf.gz
bcftools filter $sv_call_path/svaba/$sample/${sample}_svaba.vcf.gz -e 'QUAL=0' -o $sv_call_path/svaba/$sample/${sample}_svaba.qual.vcf.gz
bgzip -d $sv_call_path/svaba/$sample/${sample}_svaba.qual.vcf.gz

ml R/4.3.3
Rscript $simpleEventAnnotation \
-i $sv_call_path/svaba/$sample/${sample}_svaba.qual.vcf \
-o $sv_call_path/svaba/$sample/${sample}_svaba.qual.SEA.vcf \
-g hg38

ml purge

rm $sv_call_path/svaba/$sample/$sample.alignments.txt.gz
rm $sv_call_path/svaba/$sample/$sample.bps.txt.gz
rm $sv_call_path/svaba/$sample/$sample.discordant.txt.gz
rm $sv_call_path/svaba/$sample/$sample.log
rm $sv_call_path/svaba/$sample/$sample.contigs.bam


# cue
ref_index=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta.fai
original_path=/staging/biology/sevemthday77
change_path=/opt/cue/work
cue_vm_dir=/opt/cue
cue_sif=/opt/ohpc/Taiwania3/pkg/biology/cue/cue_v0.2.2.1.sif
ml old-module
module load libs/singularity/3.10.2
model_path=/opt/cue/data/models/cue.v2.pt

new_filename=$(echo $bam | sed "s|$original_path|$change_path|")
new_ref_index=$(echo $ref_index | sed "s|$original_path|$change_path|")

echo "bam: \"$new_filename\"" > $sv_call_path/cue/$sample/${sample}_call.data.yaml
echo "fai: \"$new_ref_index\"" >> $sv_call_path/cue/$sample/${sample}_call.data.yaml
echo "model_path: \"$model_path\"" > $sv_call_path/cue/$sample/${sample}_call.model.yaml
echo "n_cpus: $cue_thread" >> $sv_call_path/cue/$sample/${sample}_call.model.yaml

DATA_CONFIG=$(echo $sv_call_path/cue/$sample/${sample}_call.data.yaml | sed "s|$original_path|$change_path|")
MODEL_CONFIG=$(echo $sv_call_path/cue/$sample/${sample}_call.model.yaml | sed "s|$original_path|$change_path|")

singularity exec \
--bind ${original_path}:${change_path} $cue_sif \
bash -c "cd ${cue_vm_dir}; python engine/call.py --data_config ${DATA_CONFIG} --model_config ${MODEL_CONFIG}"

ml purge

rm $sv_call_path/cue/$sample/reports/predictions/*/predictions.pkl

# survivor
SURVIVOR=/opt/ohpc/Taiwania3/pkg/biology/SURVIVOR/SURVIVOR_v1.0.7/SURVIVOR

echo $sv_call_path/manta/$sample/${sample}_manta.pass.vcf  >  $sv_call_path/survivor/$sample/${sample}_vcf.files
echo $sv_call_path/delly/$sample/${sample}_delly.pass.vcf  >> $sv_call_path/survivor/$sample/${sample}_vcf.files
echo $sv_call_path/lumpy/$sample/${sample}_lumpy.vcf  >> $sv_call_path/survivor/$sample/${sample}_vcf.files
echo $sv_call_path/gridss/$sample/${sample}_gridss.pass.SEA.vcf >> $sv_call_path/survivor/$sample/${sample}_vcf.files
echo $sv_call_path/svaba/$sample/${sample}_svaba.qual.SEA.vcf  >> $sv_call_path/survivor/$sample/${sample}_vcf.files
echo $sv_call_path/cue/$sample/reports/svs.vcf >> $sv_call_path/survivor/$sample/${sample}_vcf.files

$SURVIVOR merge $sv_call_path/survivor/$sample/${sample}_vcf.files 1000 2 0 0 0 50 $sv_call_path/survivor/$sample/${sample}_sv.merged.vcf

# annotsv
PATH_raw=$PATH
module load gcc/10.5.0
export ANNOTSV=/staging/biology/sevemthday77/software/AnnotSV/AnnotSV
export PATH=/opt/ohpc/Taiwania3/pkg/biology/BEDTOOLS/BEDTOOLS_v2.29.1/bin:$PATH
export PATH=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin:$PATH
export PATH=/staging/biology/sevemthday77/software/AnnotSV/AnnotSV/bin:$PATH

AnnotSV \
-SVinputFile $sv_call_path/survivor/$sample/${sample}_sv.merged.vcf \
-outputFile $sv_call_path/annotsv/$sample/${sample}_sv.merged.annotsv.tsv \
-bcftools /opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools \
-bedtools /opt/ohpc/Taiwania3/pkg/biology/BEDTOOLS/BEDTOOLS_v2.29.1/bin/bedtools \
-genomeBuild GRCh38

ml purge
export PATH=$PATH_raw

#bam QC
fastqc=/opt/ohpc/Taiwania3/pkg/biology/FastQC/FastQC_v0.11.9/fastqc
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
samtools=/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools
Picard_path=/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.26.0

mkdir -p $bam_path/bamQC

$samtools flagstat $bam> $bam_path/bamQC/${sample}_flagstat.txt

java -jar ${Picard_path}/picard.jar CollectWgsMetrics \
I=$bam \
O=$bam_path/bamQC/${sample}_WgsMetrics.txt \
R=$ref \
INTERVALS=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38_genome_autosomal.interval_list

java -jar ${Picard_path}/picard.jar CollectGcBiasMetrics \
I=$bam \
O=$bam_path/bamQC/${sample}_GCBiasMetrics.txt \
R=$ref \
CHART=$bam_path/bamQC/${sample}_GCBiasMetrics.pdf \
S=$bam_path/bamQC/${sample}_GCBiasSummary.txt

java -jar ${Picard_path}/picard.jar CollectInsertSizeMetrics \
I=$bam \
O=$bam_path/bamQC/${sample}_InsertSizeMetrics.txt \
H=$bam_path/bamQC/${sample}_InsertSizeMetrics.pdf

java -jar ${Picard_path}/picard.jar CollectAlignmentSummaryMetrics \
I=$bam \
O=$bam_path/bamQC/${sample}_AlignmentSummaryMetrics.txt \
R=$ref

$fastqc -t $fastqc_thread -o $bam_path/bamQC $bam
