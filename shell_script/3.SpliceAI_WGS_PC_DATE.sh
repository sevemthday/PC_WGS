#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J 3.SpliceAI_WGS_PC_DATE
#SBATCH -p ngs2gpu
#SBATCH -c 12
#SBATCH --mem=180G
#SBATCH --gres=gpu:2
#SBATCH -o /home/sevemthday77/log/3.SpliceAI_WGS_PC_DATE.out.txt
#SBATCH -e /home/sevemthday77/log/3.SpliceAI_WGS_PC_DATE.err.txt
#SBATCH --mail-user=sevemthday@gmail.com
#SBATCH --mail-type=END,FAIL

# run info
date=DATE
sample=GATK_WGS_PC_${date}
base_path=/staging/biology/sevemthday77/PancreaticCancer/WGS_res

# path
spliceai_file_path=$base_path/WGS_vcf/$date/spliceai_file

# spliceai
# get vcf of new variants not in spliceai sqldb
spliceai_score_sqldb=/staging/biology/sevemthday77/software/sql_db/variants/UID2.db

get_new_vcf_variants_not_in_spliceaidb=/staging/biology/sevemthday77/software/python_/get_new_vcf_variants_not_in_spliceaidb.py
python3 $get_new_vcf_variants_not_in_spliceaidb \
--input_vcfgz $spliceai_file_path/$sample.vt.pass.pct1.vcf.gz \
--output_vcf $spliceai_file_path/$sample.vt.pass.pct1.notinspliceaidb.vcf \
--spliceai_sqldb $spliceai_score_sqldb

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
-I $spliceai_file_path/$sample.vt.pass.pct1.notinspliceaidb.vcf \
-O $spliceai_file_path/$sample.vt.pass.pct1.notinspliceaidb.spliceAI.vcf \
-R $ref \
-A $spliceai_gencode_annotation_file

bgzip $spliceai_file_path/$sample.vt.pass.pct1.notinspliceaidb.spliceAI.vcf
gatk IndexFeatureFile -I $spliceai_file_path/$sample.vt.pass.pct1.notinspliceaidb.spliceAI.vcf.gz

# update sql db
spliceai_vcfgz2sql=/staging/biology/sevemthday77/software/python_/spliceai_vcfgz2sql.py
python3 $spliceai_vcfgz2sql \
--input_spliceai_vcfgz $spliceai_file_path/$sample.vt.pass.pct1.notinspliceaidb.spliceAI.vcf.gz \
--output_spliceai_concise $spliceai_file_path/$sample.spliceAI_concise_temp \
--output_spliceai_sqldb $spliceai_score_sqldb

rm -rf $spliceai_file_path/$sample.spliceAI_concise_temp

# Get UID2 list of slivar filtered vcf
get_uid2_from_vcfgz=/staging/biology/sevemthday77/software/python_/get_uid2_from_vcfgz.py
python3 $get_uid2_from_vcfgz \
--input_vcfgz $spliceai_file_path/$sample.vt.pass.pct1.vcf.gz \
--output_UID2 $spliceai_file_path/$sample.vt.pass.pct1.vcf.UID2

# Use UID2 list to get spliceai concise txt from spliceai sql db
use_uid2_get_concise_spliceai_from_sqldb=/staging/biology/sevemthday77/software/python_/use_uid2_get_concise_spliceai_from_sqldb.py
python3 $use_uid2_get_concise_spliceai_from_sqldb \
--input_UID2 $spliceai_file_path/$sample.vt.pass.pct1.vcf.UID2 \
--output_concise_spliceai $spliceai_file_path/$sample.vt.pass.pct1.spliceai.concise.txt \
--spliceai_sqldb $spliceai_score_sqldb

ml purge
export PATH=$PATH_raw
