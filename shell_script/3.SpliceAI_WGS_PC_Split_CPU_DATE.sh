#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J 3.SpliceAI_WGS_PC_Split_CPU_DATE_FILE
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92G
#SBATCH -o /home/sevemthday77/log/3.SpliceAI_WGS_PC_Split_CPU_DATE_FILE.out.txt
#SBATCH -e /home/sevemthday77/log/3.SpliceAI_WGS_PC_Split_CPU_DATE_FILE.err.txt
#SBATCH --mail-user=sevemthday@gmail.com
#SBATCH --mail-type=END,FAIL

# run info
date=DATE
split_vcf=FILE
base_path=/staging/biology/sevemthday77/PancreaticCancer/WGS_res

# path
spliceai_file_path=$base_path/WGS_vcf/$date/spliceai_file_split

PATH_raw=$PATH
export PATH=/opt/ohpc/Taiwania3/pkg/biology/Python/Python_v2.7.18/bin:$PATH

ml old-module
ml pkg/Anaconda3
source activate /home/sevemthday77/spliceai
ml biology/GATK/4.2.0.0
ml biology/HTSLIB/1.18
ml biology/bcftools/1.13

ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
spliceai_gencode_annotation_file=/staging/biology/sevemthday77/software/spliceai/annotations/gencode.v44.basic.annotation.txt

spliceai \
-I $spliceai_file_path/$split_vcf.vcf \
-O $spliceai_file_path/$split_vcf.spliceAI.vcf \
-R $ref \
-A $spliceai_gencode_annotation_file

bgzip $spliceai_file_path/$split_vcf.spliceAI.vcf
gatk IndexFeatureFile -I $spliceai_file_path/$split_vcf.spliceAI.vcf.gz

# update sql db
spliceai_score_sqldb=/staging/biology/sevemthday77/software/sql_db/variants/UID2.db
spliceai_vcfgz2sql=/staging/biology/sevemthday77/software/python_/spliceai_vcfgz2sql.py
python3 $spliceai_vcfgz2sql \
--input_spliceai_vcfgz $spliceai_file_path/$split_vcf.spliceAI.vcf.gz \
--output_spliceai_concise $spliceai_file_path/$split_vcf.spliceAI_concise_temp \
--output_spliceai_sqldb $spliceai_score_sqldb

rm -rf $spliceai_file_path/$split_vcf.spliceAI_concise_temp

ml purge
export PATH=$PATH_raw
