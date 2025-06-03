import argparse
import csv
import pandas as pd
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        type = str, metavar = 'FILE', required = True,
        help = 'the input file with spliceai annotatio'
    )
    parser.add_argument(
        '--variant_based_output',
        type = str, metavar = 'FILE', required = True,
        help = 'the output file'
    )
    parser.add_argument(
        '--sample_based_output',
        type = str, metavar = 'FILE', required = True,
        help = 'the output file'
    )
    parser.add_argument(
        '--report_genelist',
        type = str, metavar = 'FILE', required = True,
        help = 'the gene list to report'
    )
    args = parser.parse_args()
    return args

def read_csv_as_dict(filename: str) -> dict:
    data = []
    with open(filename, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            data.append(row)
    return data

def read_report_genelist(filename: str) -> set:
    # List to store the gene names
    gene_list = []
    # Read the file
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Skip the header row and process the rest
    for line in lines[1:]:
        gene_list.append(line.strip().split('\t')[0])  # strip() removes any leading/trailing whitespace including newline
    # Print the gene list
    return set(gene_list)

def is_in_UTR_region(Impact: str) -> bool:
    UTR_impact = ["UTR5_region", "UTR3_region"]
    if "|" in Impact:
        head, tail = Impact.split("|")
        
        head_in_UTR = any(head.startswith(utr) for utr in UTR_impact)
        tail_in_UTR = any(tail.startswith(utr) for utr in UTR_impact)
        
        return head_in_UTR and tail_in_UTR
    else:
        return any(Impact.startswith(utr) for utr in UTR_impact)

def is_splicing_impact(row):
    splicing_impact = ['splice acceptor variant', 'splice donor variant']
    splicing_impact_pattern = 'coding_region_coding_variant|coding_region_splice_donor_variant|coding_region_splice_acceptor_variant'
    return row['Impact_vv'] in splicing_impact or re.search(splicing_impact_pattern, row['Impact_vv'])

def is_spliceai_variant(row):
    """
    get spliceai score report based on the given row data.

    Parameters:
    score: the column key 'spliceai_score_report' value

    Returns:
    str or int: The calculated spliceai score report.
    """
    if row['spliceai_score_report'] == '.':
        return False
    else:
        return True

def correct_ACMG(row):
    loss_of_function_genes = ['MBD4', 'POLN', 'RAD50', 'POLQ', 'DCLRE1A']
    gain_of_function_genes = ['ERBB2']

    def is_truncating_impact(row):
        truncating_impact = ['frameshift deletion', 'frameshift insertion', 'stopgain']
        return row['Impact_vv'] in truncating_impact

    def is_nmd_insensitive(row):
        return row['MANE_at_least_50nt_rule'] == 'NMD-insensitive'

    def is_truncation_less_than_ten_percent(row):
        if re.match(r'^\d', row['MANE_pct_truncated_residues']):
            return float(row['MANE_pct_truncated_residues']) < 0.1
        else:
            return False

    def lacks_domains(row):
        return not re.search('[a-zA-Z]', row['Interpro_domain']) and not re.search('[a-zA-Z]', row['unipDomain'])

    def is_pathogenic_clinvar(row):
        return (re.search('Pathogenic|pathogenic', row['CLNSIG']) and
                row['CLNREVSTAT'] in ['criteria_provided,_multiple_submitters,_no_conflicts', 'reviewed_by_expert_panel'])

    def intervar_summary_pathogenic(row):
        return row['ACMG.InterVar.summary'] in ['Pathogenic', 'Likely pathogenic']
    
    def get_intervar_correction(row):
        if row['ACMG.InterVar.summary'] == 'Pathogenic':
            return 'P_correction'
        elif row['ACMG.InterVar.summary'] == 'Likely pathogenic':
            return 'LP_correction'

    def pvs1_rule_applies(row):
        return re.search('PVS1', row['ACMG.rule.InterVar'])

    acmg_classes = row['ACMG_Classes']
    row['ACMG_correction'] = acmg_classes

    # Correct P/LP to VUS
    if acmg_classes in ['likely-pathogenic', 'pathogenic']:
        if is_truncating_impact(row):
            if is_nmd_insensitive(row) and is_truncation_less_than_ten_percent(row) and lacks_domains(row):
                row['ACMG_correction'] = 'VUS_correction'
        elif is_splicing_impact(row):
            if not is_spliceai_variant(row):
                row['ACMG_correction'] = 'VUS_correction'
        elif is_in_UTR_region(row['Impact_vv']):
            row['ACMG_correction'] = 'VUS_correction'
    
    # Correct P/LP to VUS (oncogene gene list exclude)
    if acmg_classes in ['likely-pathogenic', 'pathogenic'] and row['Symbol_vv'] in gain_of_function_genes:
        if is_truncating_impact(row):
            row['ACMG_correction'] = 'VUS_correction'
        elif is_splicing_impact(row):
            row['ACMG_correction'] = 'VUS_correction'

    # Correct VUS to P/LP (ClinVar rescue)
    if acmg_classes == 'VUS':
        if is_pathogenic_clinvar(row):
            row['ACMG_correction'] = 'P_correction'
            if row['ACMG.InterVar.summary'] == 'Pathogenic':
                row['ACMG_correction'] = 'P_correction'
            elif row['ACMG.InterVar.summary'] == 'Likely pathogenic':
                row['ACMG_correction'] = 'LP_correction'

    # Correct VUS to P/LP (InterVar rescue)
    if acmg_classes == 'VUS' and intervar_summary_pathogenic(row):
        if is_truncating_impact(row) and pvs1_rule_applies(row):
            row['ACMG_correction'] = get_intervar_correction(row)
        elif is_splicing_impact(row) and pvs1_rule_applies(row):
            if is_spliceai_variant(row):
                row['ACMG_correction'] = get_intervar_correction(row)

    # Correct VUS to P/LP (white gene list rescue)
    if acmg_classes == 'VUS' and row['Symbol_vv'] in loss_of_function_genes:
        if is_truncating_impact(row):
            if not (is_nmd_insensitive(row) and is_truncation_less_than_ten_percent(row) and lacks_domains(row)):
                row['ACMG_correction'] = 'LP_correction'
        elif is_splicing_impact(row):
            if is_spliceai_variant(row):
                row['ACMG_correction'] = 'LP_correction'

    return row


def add_filter_tag_for_report(row: dict, report_genelist) -> dict:
    """
    Process a single row of data.

    Args:
    - row (dict): A dictionary representing a row of data.
    - blacklist_UID (list): List of UID values to blacklist.
    - report_genelist (set): Dictionary containing data from a report gene list.

    Returns:
    - filtered_row (dict): The processed row with added tags.
    """

    """
    blacklist UID
    chr2-17716859-T-- (SMC6; The number of alternative reads is too low; PC0107)
    chr3-14168325---CGGCATAC (XPC; QC fail,duplicate reads; PC0097)
    chr4-186606107---CATAC (FAT1; QC fail,duplicate reads; PC0097)
    chr9-32986032-T-G (APTX; poly-A high error region; PC0073,PC0097,PC0103)
    chr22-29604105-A-G (NF2; PS1 in Deafness variant database, ClinVar:Uncertain_significance(2)|Benign(1)|Likely_benign(3); PC0007)
    """
    blacklist_UID = ["chr2-17716859-T--",
                     "chr3-14168325---CGGCATAC",
                     "chr4-186606107---CATAC",
                     "chr9-32986032-T-G",
                     "chr22-29604105-A-G"]
    
    """
    blacklist gene
    ABCB11: familial intrahepatic cholestasis [AR], Familial intrahepatic cholestasis is associated with hepatocellular carcinoma and cholangiocarcinoma [Cell.2018.cancer.predisposition.gene].
    AIP: Pituitary adenomas[AD,SMu] [CancerCell.2022.LOF.CPG.DDR.gene;DKTK.MASTER.germline142;gencc.cancer.aggregate]
    AMY2A: Amylase [pancreatic.cancer.enzyme.genelist]
    DNAJC21: cancer-prone bone marrow failure (BMF)[AR],Shwachman-Diamond syndrome[AR] -> leukemia [cancer.panel.ntuh190]
    GJB2: Deafness [Cell.2018.cancer.predisposition.gene: Keratosis-icthyosis-deafness syndrome (KID), Squamous cell carcinoma,	autosomal dominant]
    HFE: Haemochromatosis [Cell.2018.cancer.predisposition.gene: Haemochromatosis, "Hepatocellular carcinoma Cholangiocarcinoma", autosomal recessive]
    HPS1: Hermansky-Pudlak syndrome 1[AR] -> skin cancer risk [DKTK.MASTER.germline142]
    MLF1: []
    MUC5B: Idiopathic pulmonary fibrosis (IPF)[AR] -> lung cancer risk [cancer.panel.ntuh190]
    PTPRC: Immunodeficiency[AR] [cancer.gene.census.v96: TSG, HNSCC; colorectal cancer; gastric cancer; lung cancer; melanoma]
    RNF213: Moyamoya[AD,AR] [cancer.gene.census.v96: fusion, ALCL]
    SLC25A13: Citrullinaemia[AR] -> Hepatocellular carcinoma risk [Cell.2018.cancer.predisposition.gene]
    SQSTM1: osteosarcoma,Paget disease of bone 3 [cancer.panel.ntuh190;gencc.cancer.aggregate]
    TREX1: cytosolic DNA nuclease essential for regulation of cGAS-STING immune signaling [CancerCell.2022.LOF.CPG.DDR.gene;CellRep.2018.DDR.genes.276] Aberrant nuclear activity of C-terminally truncated TREX1 triggers DNA damage and senescence in mammalian cells(PMID:38824133)
    UVSSA: UV-sensitive syndrome 3[AR] [CancerCell.2022.LOF.CPG.DDR.gene;CellRep.2018.DDR.genes.276]
    VTI1A: [cancer.gene.census.v96: fusion, colorectal]
    ZNF479: [cancer.gene.census.v96: lung cancer; bladder carcinoma; prostate carcinoma]
    """
    blacklist_gene = ['ABCB11','AIP','AMY2A','DNAJC21','GJB2','HFE','HPS1','MLF1','MUC5B','PTPRC','RNF213','SLC25A13','SQSTM1','TREX1','UVSSA','VTI1A','ZNF479']

    #rescue from cancer.gene.census.v96 and gencc.cancer.aggregate
    rescue_gene = ['TRIM24','DLC1','WIF1','ACVR2A']

    tag_dict = {}

    # Add tags based on conditions
    tag_dict['tag_1'] = 'PASS' if row['QC_variant'] == 'PASS' else 'QCfail'
    tag_dict['tag_2'] = 'PASS' if (row['Panel_filtering'] not in ['cancer.gene.census.v96', 'gencc.cancer.aggregate']) or (row['Symbol_vv'] in rescue_gene) else 'Panel'
    tag_dict['tag_3'] = 'PASS' if row['TaiwanBiobank-official_Illumina1000-AF'] == '.' or float(row['TaiwanBiobank-official_Illumina1000-AF']) <= 0.01 else 'TW'
    tag_dict['tag_4'] = 'PASS' if row['Symbol_vv'] not in blacklist_gene else 'blGene'
    tag_dict['tag_5'] = 'PASS' if row['UID'] not in blacklist_UID else 'blUID'
    tag_dict['tag_6'] = 'PASS' if re.search(r'VUS|likely-pathogenic|pathogenic', row['ACMG_Classes']) else 'Benign'
    tag_dict['tag_7'] = 'PASS' if row['Symbol_vv'] in report_genelist or re.search(r'pathogenic|P_correction', row['ACMG_correction']) else 'GeneReport'
    tag_dict['tag_8'] = 'PASS' if row['rTrx'] != 'NA' else 'rTrxNA'
    tag_dict['tag_9'] = 'PASS' if row['high_recall'] != 'NA' else 'highRecallNA'
    tag_dict['tag_10'] = 'PASS' if not (is_splicing_impact(row) and not is_spliceai_variant(row)) else 'splicingVariantfail'
    tag_dict['Interpro_domain'] = row['Interpro_domain'].replace(';', '').replace('.', '')
    tag_dict['tag_11'] = 'PASS' if not (tag_dict['Interpro_domain'] == '' and row['Impact_vv'] in ['nonsynonymous SNV', 'nonframeshift insertion', 'nonframeshift deletion'] and not re.search(r'pathogenic|P_correction', row['ACMG_correction'])) else 'domain'
    tag_dict['tag_12'] = 'PASS' if not is_in_UTR_region(row['Impact_vv']) else 'UTR'
    tag_dict['tag_13'] = 'PASS' if not (row['Impact_vv'] == 'synonymous SNV' and row['spliceai_score_report'] == '.') else 'synonymousSNV'
    tag_dict['tag_14'] = 'PASS' if not (row['Impact_vv'] == 'coding_region_intron_variant' and row['spliceai_score_report'] == '.') else 'Intron'
    #tag_dict['tag_14'] = 'PASS' if row['variant_AF_avg'] != 'NA' and float(row['variant_AF_avg']) >= 0.3 else 'variant_AF_avg'
    row['filtered_variant_for_report'] = ';'.join(list(set([tag_dict[f'tag_{i}'] for i in range(1, 15)])))

    return row

def determine_zygosity(row):
    homozygous_formats = ["1|1", "1/1"]
    genotype = row['FORMAT'].split(':')[0]
    
    if genotype in homozygous_formats:
        return "Homozygous"
    else:
        return "Heterozygous"

def get_variant_report_notation(row):

    def get_t_hgvs_variant_description(row):
        if row['Zygosity'] == "Heterozygous":
            return f"c.[{row['HGVS_transcript_position_nt']}];[{row['HGVS_transcript_position']}=]"
        elif row['Zygosity'] == "Homozygous":
            return f"c.[{row['HGVS_transcript_position_nt']}];[{row['HGVS_transcript_position_nt']}]"
    
    def get_p_hgvs_variant_description(row):
        if row['HGVS_protein_short'] == '.':
            return ""
        else:
            return f"p.{row['HGVS_protein_short']}"
    
    return f"{get_t_hgvs_variant_description(row)} {get_p_hgvs_variant_description(row)}"

def convert_to_sample_based_format(variants: list) -> list:
    # get "Otherinfo12" and "Segregation" index
    column_names = list(variants[0].keys())
    index_otherinfo12 = column_names.index("Otherinfo12") # column start
    index_segregation = column_names.index("Segregation") # column end

    # Sample IDs are between "Otherinfo12" and "Segregation"
    sample_ids = column_names[index_otherinfo12 + 1: index_segregation]
    col_names_head = column_names[:index_otherinfo12 + 1]
    col_names_tail = column_names[index_segregation:]

    new_data = []
    # Split by sample ID
    for sample in sample_ids:
        for row in variants:
            if sample in row['Segregation'].split(","):
                new_row = {col: row[col] for col in col_names_head}
                new_row['FORMAT'] = row[sample]
                new_row.update({col: row[col] for col in col_names_tail})
                new_row['Sample'] = sample
                # add zygosity
                new_row['Zygosity'] = determine_zygosity(new_row)
                new_row['variant_report_notation']= get_variant_report_notation(new_row)
                new_data.append(new_row)
    return new_data

def write_csvdict_data(data: list, output: str):
    df = pd.DataFrame(data)
    df.to_csv(output, sep='\t', index=False)

def main():
    args = parse_args()
    variants = read_csv_as_dict(args.input)
    report_genelist = read_report_genelist(args.report_genelist)
    variant_based_data = []
    for row in variants:
        variant_based_data.append(add_filter_tag_for_report(correct_ACMG(row), report_genelist))
    sample_based_data = convert_to_sample_based_format(variant_based_data)
    write_csvdict_data(variant_based_data, args.variant_based_output)
    write_csvdict_data(sample_based_data, args.sample_based_output)

if __name__ == '__main__':
    main()
