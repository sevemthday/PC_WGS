########################################################################################################################
# Before manual post analysis
# 1. upload GenDiseak and download vtable
# 2. rename vtable to GATK_WGS_PC_${date}.hg38.multianno.panel.filtering_merge.vtable.tsv
# 3. upload vtable to $gene_panel_filtering_path/merge_filtered_variants/upload
########################################################################################################################

# run info
date=DATE
sample=GATK_WGS_PC_${date}
base_path=/staging/biology/sevemthday77/PancreaticCancer/WGS_res

# path
vcf_file_path=$base_path/WGS_vcf/$date/vcf_file
annotation_file_path=$base_path/WGS_vcf/$date/annotation_file
gene_panel_filtering_path=$annotation_file_path/gene_panel_filtering

#
mkdir -p $gene_panel_filtering_path/merge_filtered_variants/upload/json

# download json from variantvalidator API
download_json_from_vv_api=/staging/biology/sevemthday77/software/python_/download_json_from_vv_api.py
python3 $download_json_from_vv_api \
--UID2_input $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.UID2.clean.txt \
--json_folder $gene_panel_filtering_path/merge_filtered_variants/upload/json

# get representative transcript from vv json
#if error, add new except transcript(check in json and varsome) to appris_data.principal_MANE_except_transcript.txt
ml old-module
ml biology/Python/3.7.10
get_representative_transcript_from_vv_json=/staging/biology/sevemthday77/software/python_/get_representative_transcript_from_vv_json.py
MANEsummary=/staging/biology/sevemthday77/software/mydatabase/MANE/MANE.GRCh38.v1.1.summary.txt
python3 $get_representative_transcript_from_vv_json \
--UID2_input $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.UID2.clean.txt \
--json_folder $gene_panel_filtering_path/merge_filtered_variants/upload/json \
--output_file $gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.VV.txt \
--appris_MANE_transcript /staging/biology/sevemthday77/software/mydatabase/appris/appris_data.principal_MANE_except_transcript.txt \
--gencode_metadata_refSeq /staging/biology/sevemthday77/software/GENCODE/v44/gencode.v44.metadata.RefSeq \
--MANE_summary_txt $MANEsummary \
--refseq_gtf_file /staging/biology/sevemthday77/software/mydatabase/refseq/GCF_000001405.40_GRCh38.p14_genomic.gtf_NC.NM.exon.intron.txt

# postAnalysis
ml R/4.3.3
postAnalysisWorkflow=/staging/biology/sevemthday77/software/Rscript_/postAnalysisWorkflow.R
Rscript $postAnalysisWorkflow \
--input $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.txt \
--output $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.QC.ric.pvs1.unipDomain.txt \
--panel /staging/biology/sevemthday77/software/mydatabase/genelist/RoleInCancer/RoleInCancer.txt \
--maneInfo /staging/biology/sevemthday77/software/mydatabase/MANE/MANE.GRCh38.v1.1.summary_parseForPVS1.txt \
--unipDomainBed /staging/biology/sevemthday77/software/UCSC/hg38/uniprot/unipDomain.info.parse.bed

# merge VariantValidator and GenDiseak
mergeVariantValidatorGendiseak=/staging/biology/sevemthday77/software/Rscript_/mergeVariantValidatorGendiseak_V2.R
inputVariantValidator=$gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.VV.txt
inputGendiseak=$gene_panel_filtering_path/merge_filtered_variants/upload/$sample.hg38.multianno.panel.filtering_merge.vtable.tsv
Rscript $mergeVariantValidatorGendiseak \
--inputAnnovar $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.QC.ric.pvs1.unipDomain.txt \
--inputVariantValidator $inputVariantValidator \
--inputGendiseak $inputGendiseak \
--output $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.QC.ric.pvs1.unipDomain.vvgds.txt

# get spliceAI score report
ml old-module
ml biology/Python/3.7.10
getSpliceaiScoreDB=/staging/biology/sevemthday77/software/python_/getSpliceaiScoreDB.py
python3 $getSpliceaiScoreDB \
--input $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.QC.ric.pvs1.unipDomain.vvgds.txt \
--output $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.QC.ric.pvs1.unipDomain.vvgds.ss.txt

# Manually correct ACMG classification of variants as VUS or P/LP, and add filter tags for the report
ml old-module
ml biology/Python/3.7.10
acmgCorrectionFilter=/staging/biology/sevemthday77/software/python_/acmgCorrectionFilter.py
report_genelist=/staging/biology/sevemthday77/software/mydatabase/genelist/Report_Use_Only/NTUH.report.genelist.208.txt
python3 $acmgCorrectionFilter \
--input $gene_panel_filtering_path/merge_filtered_variants/$sample.hg38.multianno.panel.filtering_merge.QC.ric.pvs1.unipDomain.vvgds.ss.txt \
--variant_based_output $gene_panel_filtering_path/merge_filtered_variants/$sample.variant_based_final.txt \
--sample_based_output $gene_panel_filtering_path/merge_filtered_variants/$sample.sample_based_final.txt \
--report_genelist $report_genelist

# get individual report
ml old-module
ml biology/Python/3.7.10
reportOut=/staging/biology/sevemthday77/software/python_/reportOut.py
output_path=/staging/biology/sevemthday77/PancreaticCancer/WGS_res/WGS_vcf/20230621/vcf_file

python3 $reportOut \
--input $gene_panel_filtering_path/merge_filtered_variants/$sample.sample_based_final.txt \
--output $gene_panel_filtering_path/merge_filtered_variants/$sample.sample_based_final.report.txt \
--base_output_dir $gene_panel_filtering_path/merge_filtered_variants \
--sample_list_file $vcf_file_path/vcf_sample_list
