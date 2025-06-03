#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J 2.gVCF2vt_WGS_PC_DATE
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92G
#SBATCH -o /home/sevemthday77/log/2.gVCF2vt_WGS_PC_DATE.out.txt
#SBATCH -e /home/sevemthday77/log/2.gVCF2vt_WGS_PC_DATE.err.txt
#SBATCH --mail-user=sevemthday@gmail.com
#SBATCH --mail-type=END,FAIL

# run info
date=DATE
sample=GATK_WGS_PC_${date}
base_path=/staging/biology/sevemthday77/PancreaticCancer/WGS_res

# path
#fastq_path=$base_path/WGS_fastq/$date
bam_path=$base_path/WGS_bam_gvcf/$date
vcf_file_path=$base_path/WGS_vcf/$date/vcf_file
vqsr_file_path=$base_path/WGS_vcf/$date/vqsr_file
spliceai_file_path=$base_path/WGS_vcf/$date/spliceai_file

mkdir -p $vcf_file_path
mkdir -p $vqsr_file_path
mkdir -p $spliceai_file_path

# MultiQC
ml old-module
ml biology/MultiQC/1.18

cd $bam_path/bamQC
mkdir post_align_qc
cd post_align_qc
multiqc $bam_path/bamQC

ml purge

# gVCF2vt
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
GATK=/opt/ohpc/Taiwania3/pkg/biology/GATK/gatk_v4.2.0.0/gatk
slivar=/staging/biology/sevemthday77/software/slivar/slivar
gnomad=/staging/biology/sevemthday77/software/slivar/gnomad.hg38.genomes.v3.fix.zip
bcftools=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools
export PATH=/opt/ohpc/Taiwania3/pkg/biology/Python/Python_v2.7.18/bin:$PATH

ml R/4.3.3
PATH_raw=$PATH
PATH=$PATH:/opt/ohpc/Taiwania3/pkg/biology/tabix/bin

r1=/staging/biology/sevemthday77/GATK_bundle/hapmap_3.3.hg38.vcf.gz
r2=/staging/biology/sevemthday77/GATK_bundle/1000G_omni2.5.hg38.vcf.gz
r3=/staging/biology/sevemthday77/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
r4=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf
r5=/staging/biology/sevemthday77/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
r6=/staging/biology/sevemthday77/GATK_bundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
find $bam_path -name '*_raw.g.vcf.gz' > $vcf_file_path/gVCF.input.list

$GATK --java-options "-Xmx77g" CombineGVCFs \
-R $ref \
-V $vcf_file_path/gVCF.input.list \
-O $vcf_file_path/$sample.combined.g.vcf.gz

$GATK --java-options "-Xmx77g" GenotypeGVCFs \
-R $ref \
-V $vcf_file_path/$sample.combined.g.vcf.gz \
-O $vcf_file_path/$sample.joint_genotype.vcf.gz

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./" VariantRecalibrator \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 \
-tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
-R $ref \
-V $vcf_file_path/$sample.joint_genotype.vcf.gz \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${r1} \
--resource:omni,known=false,training=true,truth=true,prior=12.0 ${r2} \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${r3} \
--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${r4} \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
-mode SNP \
--max-gaussians 6 \
--rscript-file $vqsr_file_path/$sample.vqsr.snps.plots.R \
--tranches-file $vqsr_file_path/$sample.vqsr.snps.tranches \
-O $vqsr_file_path/$sample.vqsr.snps.recal

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./" ApplyVQSR \
-R $ref \
-V $vcf_file_path/$sample.joint_genotype.vcf.gz \
-mode SNP \
--recal-file $vqsr_file_path/$sample.vqsr.snps.recal \
--tranches-file $vqsr_file_path/$sample.vqsr.snps.tranches \
--truth-sensitivity-filter-level 99.7 \
--create-output-variant-index true \
-O $vcf_file_path/$sample.snps.recalibrated.vcf.gz

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./" VariantRecalibrator \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
-tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
-R $ref \
-V $vcf_file_path/$sample.snps.recalibrated.vcf.gz \
--resource:mills,known=false,training=true,truth=true,prior=12.0 ${r5} \
--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 ${r6} \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${r4} \
-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
-mode INDEL \
--max-gaussians 4 \
--rscript-file $vqsr_file_path/$sample.vqsr.indels.plots.R \
--tranches-file $vqsr_file_path/$sample.vqsr.indels.tranches \
-O $vqsr_file_path/$sample.vqsr.indels.recal

$GATK --java-options "-Xmx77G -Djava.io.tmpdir=./" ApplyVQSR \
-R $ref \
-V $vcf_file_path/$sample.snps.recalibrated.vcf.gz \
-mode INDEL \
--recal-file $vqsr_file_path/$sample.vqsr.indels.recal \
--tranches-file $vqsr_file_path/$sample.vqsr.indels.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
-O $vcf_file_path/$sample.filtered.vcf.gz

rm -f $vcf_file_path/$sample.snps.recalibrated.vcf.gz
rm -f $vcf_file_path/$sample.snps.recalibrated.vcf.gz.tbi

ml old-module
ml pkg/Anaconda3
source activate /home/sevemthday77/condapy3
vt decompose -s $vcf_file_path/$sample.filtered.vcf.gz| \
vt normalize -r $ref - | \
vt uniq - -o $vcf_file_path/$sample.vt.filtered.vcf.gz
conda deactivate
ml purge

$GATK IndexFeatureFile \
-I $vcf_file_path/$sample.vt.filtered.vcf.gz

$GATK SelectVariants \
--exclude-filtered true \
-V $vcf_file_path/$sample.vt.filtered.vcf.gz \
-O $vcf_file_path/$sample.vt.pass.vcf.gz

$slivar expr \
--vcf $vcf_file_path/$sample.vt.pass.vcf.gz \
--gnotate $gnomad \
--info "INFO.gnomad_popmax_af < 0.01" \
--pass-only \
--out-vcf $spliceai_file_path/$sample.vt.pass.pct1.vcf.gz

$bcftools index $spliceai_file_path/$sample.vt.pass.pct1.vcf.gz

ml purge
export PATH=$PATH_raw
