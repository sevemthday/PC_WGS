#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J 4.annot_WGS_PC_DATE
#SBATCH -p ngs372G
#SBATCH -c 56
#SBATCH --mem=372G
#SBATCH -o /home/sevemthday77/log/4.annot_WGS_PC_DATE.out.txt
#SBATCH -e /home/sevemthday77/log/4.annot_WGS_PC_DATE.err.txt
#SBATCH --mail-user=sevemthday@gmail.com
#SBATCH --mail-type=END,FAIL

# run info
date=DATE
sample=GATK_WGS_PC_${date}
base_path=/staging/biology/sevemthday77/PancreaticCancer/WGS_res

# path
vcf_file_path=$base_path/WGS_vcf/$date/vcf_file
spliceai_file_path=$base_path/WGS_vcf/$date/spliceai_file
annotation_file_path=$base_path/WGS_vcf/$date/annotation_file
intervar_file_path=$base_path/WGS_vcf/$date/intervar_file
gene_panel_filtering_path=$annotation_file_path/gene_panel_filtering

#
mkdir -p $annotation_file_path
mkdir -p $intervar_file_path

#
ref=/staging/biology/sevemthday77/GATK_bundle/Homo_sapiens_assembly38.fasta
table_annovar=/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl
ANNOVAR_humandb=/staging/biology/sevemthday77/software/mydatabase/ANNOVAR/humandb/
bcftools=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools
InterVar_path=/staging/biology/sevemthday77/software/InterVar/InterVar-2.2.1/Intervar.py
intervardb=/staging/biology/sevemthday77/software/InterVar/InterVar-2.2.1/intervardb
ANNOVAR_path=/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819
export PATH=/opt/ohpc/Taiwania3/pkg/biology/Python/Python_v2.7.18/bin:$PATH

$bcftools query -l $vcf_file_path/$sample.vt.filtered.vcf.gz > $vcf_file_path/vcf_sample_list

$table_annovar \
$vcf_file_path/$sample.vt.filtered.vcf.gz \
$ANNOVAR_humandb \
-buildver hg38 \
-out $annotation_file_path/$sample \
-remove \
-protocol \
refGene,knownGene,ensGene,clinvar_20221231,clinvar_CLNSIGCONF_20221231,TWB1492_AF,TaiwanBiobank-official,gnomad30_genome,gnomad211_genome,gnomad211_exome,avsnp150,dbnsfp42a,cosmic94_coding,cosmic94_noncoding,icgc28,dbscsnv11 \
-operation gx,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f \
-arg '-splicing 10',,,,,,,,,,,,,,, \
-nastring . \
-vcfinput \
-polish \
--maxgenethread 24 \
--thread 20

# "After ANNOVAR annotation, you must immediately perform sample segregation, otherwise there will be errors when outputting the VCF later."
ml R/4.3.3

sampleSegregation=/staging/biology/sevemthday77/software/Rscript_/sampleSegregation.R
Rscript $sampleSegregation \
-i $annotation_file_path/$sample.hg38_multianno.txt \
-o $annotation_file_path/$sample.hg38_multianno.segregation.txt \
--sampleID $vcf_file_path/vcf_sample_list

get_uid2_in_annovartxt=/staging/biology/sevemthday77/software/python_/get_uid2_in_annovartxt.py
python3 $get_uid2_in_annovartxt \
--input $annotation_file_path/$sample.hg38_multianno.segregation.txt \
--output $annotation_file_path/$sample.hg38_multianno.segregation.uid2.txt

CLNSIGCONFmajor=/staging/biology/sevemthday77/software/Rscript_/CLNSIGCONFmajor.R
Rscript $CLNSIGCONFmajor \
-i $annotation_file_path/$sample.hg38_multianno.segregation.uid2.txt \
-o $annotation_file_path/$sample.hg38_multianno.segregation.uid2.CLNSIGCONFmajor.txt

# intervar
python $InterVar_path \
-i $annotation_file_path/$sample.avinput \
-o $intervar_file_path/$sample.intervar \
-b hg38 \
-t $intervardb \
--table_annovar=$ANNOVAR_path/table_annovar.pl \
--convert2annovar=$ANNOVAR_path/convert2annovar.pl \
--annotate_variation=$ANNOVAR_path/annotate_variation.pl \
-d $ANNOVAR_humandb

InterVarAggregate=/staging/biology/sevemthday77/software/Rscript_/InterVarAggregate.R
Rscript $InterVarAggregate \
-i $intervar_file_path/$sample.intervar.hg38_multianno.txt.intervar \
-o $intervar_file_path/$sample.intervar.hg38_multianno.txt.intervar.aggregate

mergeTableViaUID=/staging/biology/sevemthday77/software/Rscript_/mergeTableViaUID.R
Rscript $mergeTableViaUID \
-i $annotation_file_path/$sample.hg38_multianno.segregation.uid2.CLNSIGCONFmajor.txt \
-I $intervar_file_path/$sample.intervar.hg38_multianno.txt.intervar.aggregate \
-o $annotation_file_path/$sample.hg38_multianno.segregation.uid2.CLNSIGCONFmajor.mergeIntervar.txt

#left joint
mergeTwoTables=/staging/biology/sevemthday77/software/python_/mergeTwoTables.py
python3 $mergeTwoTables \
--input_file_1 $annotation_file_path/$sample.hg38_multianno.segregation.uid2.CLNSIGCONFmajor.mergeIntervar.txt \
--input_file_2 $spliceai_file_path/$sample.vt.pass.pct1.spliceai.concise.txt \
--output_file $annotation_file_path/$sample.hg38_multianno.final.txt \
--temp_sqldb $annotation_file_path/$sample.merge.temp.db \
--key_column UID2


#
panelDir=/staging/biology/sevemthday77/software/mydatabase/genelist/PC_WGS

#
mkdir -p $gene_panel_filtering_path
mkdir -p $gene_panel_filtering_path/merge_filtered_variants
mkdir -p $gene_panel_filtering_path/merge_filtered_variants/upload

#
ml R/4.3.3
annovarFilterCombinePanel=/staging/biology/sevemthday77/software/Rscript_/annovarFilterCombinePanel.R

Rscript $annovarFilterCombinePanel \
-i $annotation_file_path/$sample.hg38_multianno.final.txt \
-o $gene_panel_filtering_path/$sample.hg38_multianno \
--gnomad TRUE \
--twbiobank TRUE \
--gnomadCutoff 0.05 \
--twbiobankCutoff 0.05 \
--regionRefGene exonic:splicing \
--remvSynonymousSNV TRUE \
--keepPASS TRUE \
--panelDir $panelDir \
--clinvarPLP TRUE

Rscript $annovarFilterCombinePanel \
-i $annotation_file_path/$sample.hg38_multianno.final.txt \
-o $gene_panel_filtering_path/$sample.pct1.spliceAI.hg38_multianno \
--gnomad TRUE \
--twbiobank TRUE \
--gnomadCutoff 0.01 \
--twbiobankCutoff 0.01 \
--keepPASS TRUE \
--panelDir $panelDir \
--spliceai TRUE \
--spliceaiCutoff recommended

# merge variants and generate VCF to upload GenDiseak
mergeTable=/staging/biology/sevemthday77/software/Rscript_/mergeTable.R

Rscript $mergeTable \
-I *_panel.filtering.txt$ \
--path $gene_panel_filtering_path \
-O $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.txt \
--dedup \
--col_char

# VCF
annovarSubsetForVCF=/staging/biology/sevemthday77/software/Rscript_/annovarSubsetForVCF.R
Rscript $annovarSubsetForVCF \
-I $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.txt \
-O $gene_panel_filtering_path/merge_filtered_variants/upload/temp.txt

ml old-module
ml biology/bcftools/1.13
bcftools view -h $vcf_file_path/$sample.vt.filtered.vcf.gz > $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.vcf
cat $gene_panel_filtering_path/merge_filtered_variants/upload/temp.txt >> $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.vcf
rm -rf $gene_panel_filtering_path/merge_filtered_variants/upload/temp.txt

# Pseudo VCF
annovarSubsetForPseudoVCF=/staging/biology/sevemthday77/software/Rscript_/annovarSubsetForPseudoVCF.R
Rscript $annovarSubsetForPseudoVCF \
-I $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.txt \
-O $gene_panel_filtering_path/merge_filtered_variants/upload/temp.txt

ml old-module
ml biology/bcftools/1.13
bcftools view -h $vcf_file_path/$sample.vt.filtered.vcf.gz|grep -P '^##' > $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.pseudo.vcf
cat $gene_panel_filtering_path/merge_filtered_variants/upload/temp.txt >> $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.pseudo.vcf
rm -rf $gene_panel_filtering_path/merge_filtered_variants/upload/temp.txt

ml purge

# annovar subset for vv
ml R/4.3.3
annovarSubsetForVV=/staging/biology/sevemthday77/software/Rscript_/annovarSubsetForVV.R
Rscript $annovarSubsetForVV \
-I $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.txt \
-O $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.UID2.txt

grep -v '\*' $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.UID2.txt > $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.UID2.clean.txt
