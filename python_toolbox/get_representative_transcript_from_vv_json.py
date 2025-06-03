import argparse
import csv
import json
import pandas as pd
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--UID2_input',
        type = str, metavar = 'FILE', required = True,
        help = 'VCF-like variants: chr-pos-ref-alt'
    )
    parser.add_argument(
        '--json_folder',
        type = str, metavar = 'FOLDER', required = True,
        help = 'the path of json files'
    )
    parser.add_argument(
        '--output_file',
        type = str, metavar = 'FILE', required = True,
        help = 'the path of output file'
    )
    parser.add_argument(
        '--appris_MANE_transcript',
        type = str, metavar = 'FILE', required = True,
        help = 'the appris and MANE transcript with a few exception'
    )
    parser.add_argument(
        '--gencode_metadata_refSeq',
        type = str, metavar = 'FILE', required = True,
        help = 'the gencode metadata with ensembl transcript and its paired refseq'
    )
    parser.add_argument(
        '--MANE_summary_txt',
        type = str, metavar = 'FILE', required = True,
        help = 'the gencode metadata with ensembl transcript and its paired refseq'
    )
    parser.add_argument(
        '--refseq_gtf_file',
        type = str, metavar = 'FILE', required = True,
        help = 'the parsed refseq gtf with exon/intron region and exon_number information'
    )
    args = parser.parse_args()
    return args

# Read the contents of the file and obtain variants
def read_file_and_get_variants(UID2_input: str) -> list:
    """
    read UID2 variant file 

    Args:
        UID2_input: vcf-like variants with structure "chrom-pos-ref-alt" in the file

    Returns:
        list: The list contain the UID2 variants
    """
    with open(UID2_input, "r") as f:
        variants = [line.strip() for line in f if line.strip()]
    return variants

def parse_appris_transcript(transcript_file: str) -> set:
    """
    read the transcripts from appris, MANE and except_transcript

    Args:
        transcript_file: transcript list per line

    Returns:
        set: a set of transcripts
    """
    with open(transcript_file, 'r') as f:
        transcript_list = []
        for line in f:
            transcript = line.strip().split(".")[0]
            transcript_list.append(transcript)
    return(set(transcript_list))

def parse_gencode_nm_to_enst(file_path: str) -> dict:
    '''
    read the gencode metadata RefSeq which contain the ensembl transcript id and its associated refseq

    Args:
        file_path: paired transcripts (ensembl and refseq) list per line

    Returns:
        dict with key refseq transcript ID and with value ensembl transcript id
    '''
    nm_to_enst = {}
    with open(file_path, 'r') as f:
        for line in f:
            line_split = re.split(r'\s+', line.strip())
            if len(line_split) >= 2:
                enst, nm = line_split[:2]
                #enst = enst.split(".")[0]
                nm = nm.split(".")[0]
                nm_to_enst.setdefault(nm, set()).add(enst)
    # retun: a dic with key nm and value enst_set
    return {nm: list(enst_set) for nm, enst_set in nm_to_enst.items()}

def get_matching_transcripts(selected_Trx: dict, Trx_all: list) -> list:
    """
    Retrieves matching transcripts based on HGVS notation (selected_Trx) from a list of transcripts (Trx_all).

    Args:
        selected_Trx (dict): A dictionary representing the selected transcript with HGVS notation.
        Trx_all (list): A list of dictionaries representing all transcripts.

    Returns:
        list: A list containing combined matching transcripts in NM and ENST formats.
    """
    # Initialize lists to store matching transcripts
    matching_transcripts_nm = []
    matching_transcripts_enst = []

    # Iterate through each transcript in the list
    for row in Trx_all:
        # Check if the transcript has a HGVS notation
        if row["T_HGVS"] is not None:
            # Split the transcript HGVS to get the accession number
            split_transcript_hgvs = row["T_HGVS"].split(":")
            # Check if the HGVS notation matches the selected transcript's HGVS notation
            if split_transcript_hgvs[1] == selected_Trx["T_HGVS"].split(":")[1]:
                # Check if the transcript starts with NM_ (RefSeq) or ENST (Ensembl)
                if row["T_HGVS"].startswith("NM_"):
                    matching_transcripts_nm.append(row["rTrx"])
                elif row["T_HGVS"].startswith("ENST"):
                    matching_transcripts_enst.append(row["rTrx"])

    # Combine matching transcripts into strings or set defaults if no matches found
    combined_transcripts_nm = ",".join(matching_transcripts_nm) if matching_transcripts_nm else "."
    combined_transcripts_enst = ",".join(matching_transcripts_enst) if matching_transcripts_enst else "."

    return [combined_transcripts_nm, combined_transcripts_enst]

def find_paired_transcripts(nm_number: str, nm_to_enst: dict) -> list:
    '''
    Function to find paired ENST numbers for a given NM number
    '''
    if nm_number in nm_to_enst:
        return nm_to_enst[nm_number]
    else:
        return ['.']

def extract_csv_data(csv_file, key_column, value_columns, condition_column=None, condition_value=None):
    """
    Extracts data from a CSV file and organizes it into a dictionary.

    Args:
        csv_file (str): The path to the CSV file.
        key_column (str): The column name to use as keys in the result dictionary.
        value_columns (list): A list of column names whose values will be extracted and stored under each key.
        condition_column (str, optional): The column name for the condition check. Defaults to None.
        condition_value (str, optional): The value to check in the condition_column. Defaults to None.

    Returns:
        dict: A dictionary containing extracted data organized by keys from key_column.
    """
    result_dict = {}
    with open(csv_file, newline='') as f:
        csv_reader = csv.DictReader(f, delimiter='\t')
        for row in csv_reader:
            # Check if a condition column and value are provided, and if the condition matches
            if condition_column is None or row[condition_column] == condition_value:
                key_value = row[key_column]
                # Initialize an empty list if the key is not yet in the result dictionary
                if key_value not in result_dict:
                    result_dict[key_value] = []
                # Extract values from specified columns and extend the values list for the key
                values = [row[col] for col in value_columns]
                result_dict[key_value].extend(values)
    return result_dict

def parse_mane_summary_txt(file_path: str, condition_value: str) -> dict:
    """
    Parses a MANE summary TXT file and extracts relevant data based on a condition.

    Args:
        file_path (str): The path to the MANE summary TXT file.
        condition_value (str): The value to check in the 'MANE_status' column.

    Returns:
        dict: A dictionary containing extracted data organized by 'symbol' with 'RefSeq_nuc' and 'Ensembl_nuc' values.
    """
    # Extract data using the extract_csv_data function based on the MANE status condition
    result_dict = extract_csv_data(file_path, 'symbol', ['RefSeq_nuc', 'Ensembl_nuc'], 'MANE_status', condition_value)
    return result_dict

def get_symbol_paired_transcripts(symbol: str, symbol_to_trx: dict) -> list:
    """
    Retrieves paired transcripts ['RefSeq_nuc', 'Ensembl_nuc'] for a given symbol from a dictionary.

    Args:
        symbol (str): The symbol for which paired transcripts are requested.
        symbol_to_trx (dict): A dictionary mapping symbols to lists of transcripts.

    Returns:
        list: A list containing paired transcripts for the symbol, or ['.', '.'] if the symbol is not found.
    """
    # Check if the symbol exists in the symbol_to_trx dictionary
    if symbol in symbol_to_trx:
        # Return the paired transcripts if the symbol is found
        return symbol_to_trx[symbol]
    else:
        # Return ['.', '.'] if the symbol is not found
        return ['.', '.']

def analyze_intron_variant(variant_str):
    # Step 1: Check if there is an underscore (_) after the "c."
    variant_str = variant_str.split(":", 1)[1]
    if "_" in variant_str:
        # Step 2: Split the string into two parts based on the underscore (_)
        part1, part2 = variant_str.split("_", 1)
        result_part1 = analyze_c_part1(part1)
        result_part2 = analyze_c_part2(part2)
        if result_part1 == result_part2:
            return result_part1
        else:
            return f"{result_part1}|{result_part2}"
    else:
        return analyze_c_part1(variant_str)

def analyze_c_part1(c_part_str):
    # Step 3: Check for specific patterns in the c. part
    if re.match(r"c.-\d+[+][12](?!\d)", c_part_str):
        return "UTR5_region_splice_donor_variant"
    elif re.match(r"c.-\d+[-][12](?!\d)", c_part_str):
        return "UTR5_region_splice_acceptor_variant"
    elif re.match(r"c.-\d+[+-]", c_part_str):
        return "UTR5_region_intron_variant"
    elif re.match(r"c.-\d+", c_part_str):
        return "UTR5_region_UTR5_variant"
    elif re.match(r"c.\*\d+[+][12](?!\d)", c_part_str):
        return "UTR3_region_splice_donor_variant"
    elif re.match(r"c.\*\d+[-][12](?!\d)", c_part_str):
        return "UTR3_region_splice_acceptor_variant"
    elif re.match(r"c.\*\d+[+-]", c_part_str):
        return "UTR3_region_intron_variant"
    elif re.match(r"c.\*\d+", c_part_str):
        return "UTR3_region_UTR3_variant"
    elif re.match(r"c.\d+[+][12](?!\d)", c_part_str):
        return "coding_region_splice_donor_variant"
    elif re.match(r"c.\d+[-][12](?!\d)", c_part_str):
        return "coding_region_splice_acceptor_variant"
    elif re.match(r"c.\d+[+-]", c_part_str):
        return "coding_region_intron_variant"
    elif re.match(r"c.\d+", c_part_str):
        return "coding_region_coding_variant"

def analyze_c_part2(c_part_str):
    # Step 3: Check for specific patterns in the c. part
    if re.match(r"-\d+[+][12](?!\d)", c_part_str):
        return "UTR5_region_splice_donor_variant"
    elif re.match(r"-\d+[-][12](?!\d)", c_part_str):
        return "UTR5_region_splice_acceptor_variant"
    elif re.match(r"-\d+[+-]", c_part_str):
        return "UTR5_region_intron_variant"
    elif re.match(r"-\d+", c_part_str):
        return "UTR5_region_UTR5_variant"
    elif re.match(r"\*\d+[+][12](?!\d)", c_part_str):
        return "UTR3_region_splice_donor_variant"
    elif re.match(r"\*\d+[-][12](?!\d)", c_part_str):
        return "UTR3_region_splice_acceptor_variant"
    elif re.match(r"\*\d+[+-]", c_part_str):
        return "UTR3_region_intron_variant"
    elif re.match(r"\*\d+", c_part_str):
        return "UTR3_region_UTR3_variant"
    elif re.match(r"\d+[+][12](?!\d)", c_part_str):
        return "coding_region_splice_donor_variant"
    elif re.match(r"\d+[-][12](?!\d)", c_part_str):
        return "coding_region_splice_acceptor_variant"
    elif re.match(r"\d+[+-]", c_part_str):
        return "coding_region_intron_variant"
    elif re.match(r"\d+", c_part_str):
        return "coding_region_coding_variant"

def classify_hgsv_variant(p_vcf, t_hgvs, p_hgvs_tlc):
    '''
    Categorize HGVS variants based on their coding impact.
    p_vcf: pseudo vcf
    t_hgvs: HGSV transcript
    p_hgvs_tlc: HGSV protein 3-letter
    splicing donor (+1 and +2) and acceptor (-1 and -2)
    '''
    # parse p_vcf
    ref, alt = p_vcf.split('-')[-2:]  # Extract reference and alternative nucleotides
    ref_len = len(ref)
    alt_len = len(alt)
    #
    try:
        if re.match(r'NP_\d+\.\d+:p\.\([A-Za-z]+\d+[A-Za-z]+fsTer\d+\)', p_hgvs_tlc):
            if ref_len > alt_len:
                return "frameshift deletion"
            elif ref_len < alt_len:
                return "frameshift insertion"
            else:
                return "unknown"
        elif re.match(r'NP_\d+\.\d+:p\.\([A-Za-z]+\d+Ter\)', p_hgvs_tlc):
            return "stopgain"
        elif re.match(r'NP_\d+\.\d+:p\.\(Ter\d+[A-Za-z]+extTer\d+\)', p_hgvs_tlc):
            return "stoploss"
        elif re.match(r'NP_\d+\.\d+:p\.\([A-Za-z0-9_]+dup\)', p_hgvs_tlc):
            return "nonframeshift insertion"
        elif re.match(r'NP_\d+\.\d+:p\.\([A-Za-z0-9_]+ins[A-Za-z]+\)', p_hgvs_tlc):
            return "nonframeshift insertion"
        elif re.match(r'NP_\d+\.\d+:p\.\([A-Za-z0-9_]+del\)', p_hgvs_tlc):
            return "nonframeshift deletion"
        elif re.match(r'NP_\d+\.\d+:p\.\([A-Za-z]+\d+[A-Za-z]+\)', p_hgvs_tlc):
            return "nonsynonymous SNV"
        elif re.match(r'NP_\d+\.\d+:p\.\([A-Za-z]+\d+=\)', p_hgvs_tlc):
            return "synonymous SNV"
        elif re.match(r'NP_\d+\.\d+:p\.\(Met\d+\?\)', p_hgvs_tlc):
            return "startloss"
        elif re.match(r'NP_\d+\.\d+:p\.\?', p_hgvs_tlc):
            if re.match(r"NM_\d+\.\d+:c\.\d+[+][12][A-Z]>", t_hgvs):
                return "splice donor variant"
            elif re.match(r"NM_\d+\.\d+:c\.\d+[-][12][A-Z]>", t_hgvs):
                return "splice acceptor variant"
            else:
                return analyze_intron_variant(t_hgvs)
        else:
            return "unknown"
    except TypeError:
        return "unknown"

def fetch_json_transcript_to_table(vv_json_data: dict, variant_description: str, appris_Trx: set) -> list:
    # Extract data from hgvs_t_and_p and convert it into a table
    # Trx: save wanted data. Trx_all: to check data only
    Impact_precedence = {
    'frameshift insertion': 1,
    'frameshift deletion': 2,
    'stopgain': 3,
    "splice donor variant": 4,
    "splice acceptor variant": 5,
    'startloss': 6,
    'stoploss': 7,
    'nonframeshift insertion': 8,
    'nonframeshift deletion': 9,
    'nonsynonymous SNV': 10,
    'synonymous SNV': 11,
    'intron variant': 12,
    'unknown': 13
    }
    Trx = []
    Trx_all = []
    for transcript, variant_data in vv_json_data[variant_description][variant_description]["hgvs_t_and_p"].items():
        gene_info = variant_data.get("gene_info", {"hgnc_id": "", "symbol": ""})
        select_status = variant_data.get("select_status", {})
        row = {
            "UID2": variant_description,
            "rTrx": transcript,
            "T_HGVS": variant_data.get("t_hgvs", ""),
            "P_HGVS_SLC": variant_data.get("p_hgvs_slc", ""),
            "P_HGVS_TLC": variant_data.get("p_hgvs_tlc", ""),
            "Mane_Select": select_status.get("mane_select", False) if select_status is not None else False,
            "Symbol_vv": gene_info["symbol"],
            "Impact_vv": classify_hgsv_variant(variant_description, variant_data.get("t_hgvs", ""), variant_data.get("p_hgvs_tlc", ""))
        }
        transcript_split_head = transcript.split(".")[0]
        if transcript_split_head in appris_Trx:
            Trx.append(row)
        Trx_all.append(row)
    print(f'variant: {variant_description}')
    print(f'all representative transcript info: {Trx}')
    sorted_Trx = sorted(Trx, key=lambda x: (Impact_precedence.get(x['Impact_vv'], 13), not x['Mane_Select']))
    return [sorted_Trx, Trx_all]

def get_representative_transcript(variants, json_file_path, appris_Trx, output_file, nm_to_enst, symbol_to_MANEtrx, symbol_to_MANEtrx_clinical, refseq_gtf_dict):
    '''
    1. get representative transcript according to appris, MANE and a few excepted transcript.
    2. The excepted transcripts are used because those transcripts, which used in VariantValidator are not in appris and MANE database.
    3. The excepted transcripts are selected based on varsome(major), VariantValidator(major) and ANNOVAR(partial).
    4. remove VariantValidator HGSV: 'intergenic' and 'NR_*'
    '''
    # json parse
    # rTrx: representative transcript
    # variants: a set of variants from UID2; variant_data: transcript-based hgsv info;
    rTrx = []
    for variant_description in variants:
        with open(f'{json_file_path}/{variant_description}.json', 'r') as json_file:
            vv_json_data = json.load(json_file)
        # Extract data from hgvs_t_and_p and convert it into a table
        # Trx: save wanted data. Trx_all: to check data only
        sorted_Trx, Trx_all = fetch_json_transcript_to_table(vv_json_data, variant_description, appris_Trx)
        # Function to get the highest impact transcript considering MANE
        check_intergenic = len(Trx_all) == 1 and Trx_all[0]['rTrx'] == 'intergenic'
        # Assuming Trx_all is a list of dictionaries
        check_all_nrTrx = all(item['rTrx'].startswith('NR_') or item['rTrx'].startswith('ENST') for item in Trx_all)
        if check_intergenic:
            pass  # skip appending
        elif check_all_nrTrx:
            pass
        else:
            selected_Trx = sorted_Trx[0]
            # get matching transcripts nm and enst if their HGVS_transcript_nucleotide_change are the same with rTrx
            combined_transcripts_nm, combined_transcripts_enst = get_matching_transcripts(selected_Trx, Trx_all)
            selected_Trx["matching_rTrx_nm"] = combined_transcripts_nm
            selected_Trx["matching_rTrx_enst"] = combined_transcripts_enst
            # get paired ensembl transcript according to GENCODE meta data
            paired_ensembl_Trx = find_paired_transcripts(selected_Trx["rTrx"].split(".")[0], nm_to_enst)
            selected_Trx["rTrx_paired_enst"] = ",".join(paired_ensembl_Trx)
            # get paired MANE refseq and ensembl transcript according to symbol of MANE.GRCh38.v1.1.summary.txt
            selected_Trx["t_mane_select_refseq"] = get_symbol_paired_transcripts(selected_Trx["Symbol_vv"], symbol_to_MANEtrx)[0]
            selected_Trx["t_mane_plus_clinical_refseq"] = get_symbol_paired_transcripts(selected_Trx["Symbol_vv"], symbol_to_MANEtrx_clinical)[0]
            selected_Trx["t_mane_select_ensembl"] = get_symbol_paired_transcripts(selected_Trx["Symbol_vv"], symbol_to_MANEtrx)[1]
            selected_Trx["t_mane_plus_clinical_ensembl"] = get_symbol_paired_transcripts(selected_Trx["Symbol_vv"], symbol_to_MANEtrx_clinical)[1]
            # parse hgvs transcript and protein into small elements
            selected_Trx = parse_hgvs_transcript_protein(selected_Trx)
            # get exon and intron number of variant
            selected_Trx = get_exon_intron_number_by_check_overlap_gtf(selected_Trx, refseq_gtf_dict)
            rTrx.append(selected_Trx)
    df = pd.DataFrame(rTrx)
    df.to_csv(output_file, sep='\t', index=False)

def get_t_hgvs_nucleotide_change(variant: str) -> str:
    # Step 1: Split based on "c." and select the second part
    nucleotide_change = variant.split("c.")[1]
    return nucleotide_change

def get_t_hgvs_position(variant: str) -> str:
    '''
    Other programming approaches:
    def extract_variant_position(variant):
    match = re.search(r'c\.(\d+[+-]?\d*[_+]?\d*|\d+_\d+)(\+\d+)?\w*', variant)
    if match:
        return match.group(1) + (match.group(2) if match.group(2) else '')
    else:
        return None
    '''
    # Step 1: Split based on "c." and select the second part
    nucleotide_change = variant.split("c.")[1]
    # Step 2: Split each result from Step 1 based on letters [a-zA-Z] and select the first part
    extracted_position = re.split('[a-zA-Z]+', nucleotide_change)[0]
    return extracted_position

def get_p_hgvs_slc_short(variant: str) -> str:
    protein_change = variant.split(":p.")[1]
    if protein_change == "?":
        return "."
    else:
        protein_change = protein_change.replace("(", "").replace(")", "")
        return protein_change

def get_p_hgvs_tlc_position(variant: str) -> str:
    # Split by ":p."
    protein_change = variant.split(":p.")[1]
    if protein_change == "?":
        return "."
    else:   
        numbers = re.findall(r'\d+', protein_change)
        if "_" in protein_change:
            return f"{numbers[0]}_{numbers[1]}"
        else:
            return numbers[0]

def get_p_hgvs_tlc_ref(variant: str, variant_impact: str) -> str:
    # Split by ":p."
    protein_change = variant.split(":p.")[1]
    if protein_change == "?":
        return "."
    else:
        protein_change = protein_change.replace("(", "").replace(")", "")
        impact_set = {"frameshift deletion", "frameshift insertion", "nonsynonymous SNV", "stopgain", "synonymous SNV"}
        impact_set_escape = {"nonframeshift deletion", "nonframeshift insertion", "startloss", "stoploss"}
        if variant_impact in impact_set:
            return protein_change[:3]
        elif variant_impact in impact_set_escape:
            return "."

def get_p_hgvs_tlc_alt(variant: str, variant_impact: str) -> str:
    # Split by ":p."
    protein_change = variant.split(":p.")[1]
    if protein_change == "?":
        return "."
    else:
        protein_change = protein_change.replace("(", "").replace(")", "")
        impact_set_escape = {"nonframeshift deletion", "nonframeshift insertion", "startloss", "stopgain", "stoploss", "synonymous SNV"}
        if variant_impact in {"frameshift deletion", "frameshift insertion"}:
            letters_before_fsTer = protein_change.split("fsTer")[0]
            return letters_before_fsTer[-3:]
        elif variant_impact == "nonsynonymous SNV":
            return protein_change[-3:]
        elif variant_impact in impact_set_escape:
            return "."

def get_p_hgvs_tlc_position_ter(variant: str, variant_impact: str) -> str:
    # Split by ":p."
    protein_change = variant.split(":p.")[1]
    if protein_change == "?":
        return "."
    else:
        protein_change = protein_change.replace("(", "").replace(")", "")
        if "Ter" in protein_change:
            position_ter = protein_change.split("Ter")[1]
            if position_ter:
                return position_ter
            else:
                return "."
        else:
            return "."

def parse_hgvs_transcript_protein(row: dict) -> dict:
    # HGVS_transcript_position_nt: string, HGVS_transcript_position: int
    # HGVS_protein_short: string, HGVS_protein_position: int, HGVS_protein_ref: tlc, HGVS_protein_alt: tlc, HGVS_Protein_position_ter: int
    row["HGVS_transcript_position_nt"] = get_t_hgvs_nucleotide_change(row["T_HGVS"])
    row["HGVS_transcript_position"] = get_t_hgvs_position(row["T_HGVS"])
    row["HGVS_protein_short"] = get_p_hgvs_slc_short(row["P_HGVS_SLC"])
    row["HGVS_protein_position"] = get_p_hgvs_tlc_position(row["P_HGVS_TLC"])
    row["HGVS_protein_ref"] = get_p_hgvs_tlc_ref(row["P_HGVS_TLC"], row["Impact_vv"])
    row["HGVS_protein_alt"] = get_p_hgvs_tlc_alt(row["P_HGVS_TLC"], row["Impact_vv"])
    row["HGVS_Protein_position_ter"] =get_p_hgvs_tlc_position_ter(row["P_HGVS_TLC"], row["Impact_vv"])
    return row

def parse_gtf_to_transcript_dict(csv_file: str, columns_to_keep: list) -> dict:
    '''
    csv columns: seqnames, start, end, type, transcript_id, exon_number
    key: transcript_id
    '''
    data = []
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            data.append(row)
    transcript_dict = {}
    for row in data:
        transcript_id = row['transcript_id']
        transcript_id = transcript_id.split(".")[0]
        row_dict = {key: row[key] for key in row if key in columns_to_keep}
        if transcript_id in transcript_dict:
            transcript_dict[transcript_id].append(row_dict)
        else:
            transcript_dict[transcript_id] = [row_dict]
    return transcript_dict

def get_variant_start_end(UID2: str)-> tuple:
    parts = UID2.split('-')
    chromosome = parts[0]
    if len(parts[2]) > len(parts[3]):
        start = int(parts[1]) + 1
        end = start + len(parts[2]) - 2  # Calculate end position based on length of ALT allele
    else:
        start = end = int(parts[1])
    #-> tuple[str, int, int]
    return chromosome, start, end

def get_exon_intron_number_by_check_overlap_gtf(row, trx_dict):
    '''
    Annotates a variant with exon and intron numbers based on overlapping regions in a GTF file.
    Args:
        row (dict): A dictionary representing a row of data containing the variant and transcript information.
        transcript_dict (dict): A dictionary containing transcript information keyed by transcript IDs.
    Returns:
        dict: The input row dictionary annotated with exon and intron numbers.
    trx: transcript
    To-do-list: variant (chr17-82009034-C-CAGCCCGTGGACCGGG) NM_024083.4:c.951_965dup	NP_076988.1:p.(D319_V323dup) but report position in rTrx_intron_no:7
    '''
     # Extract necessary information from the row
    variant = row['UID2']
    trx_id = row['rTrx'].split(".")[0]
    variant_chrom, variant_start, variant_end = get_variant_start_end(variant)

    # Check if trx_id exists in trx_dict
    if trx_id not in trx_dict:
        # Set default values for exon and intron numbers
        row["rTrx_exon_no"] = '.'
        row["rTrx_intron_no"] = '.'
        return row

    # Get the GTF information for the transcript
    trx_gtf = trx_dict[trx_id]

    # Find overlapping regions
    overlap_regions = []
    for region in trx_gtf:
        region_chrom = region['seqnames']
        region_start = int(region['start'])
        region_end = int(region['end'])
        
        # Check for overlap
        if variant_chrom == region_chrom and variant_start <= region_end and variant_end >= region_start:
            overlap_regions.append({"type": region['type'], "number": region['exon_number']})
    
    # Initialize exon and intron numbers
    row["rTrx_exon_no"] = '.'
    row["rTrx_intron_no"] = '.'

    # Update exon and intron numbers based on overlap regions
    for overlap_region in overlap_regions:
        if overlap_region["type"] == "exon":
            row["rTrx_exon_no"] = overlap_region["number"]
        elif overlap_region["type"] == "intron":
            row["rTrx_intron_no"] = overlap_region["number"]

    return row

if __name__ == '__main__':
    args = parse_args()
    UID2_input = args.UID2_input
    json_file_path = args.json_folder
    output_file = args.output_file
    appris_MANE_transcript_file = args.appris_MANE_transcript
    gencode_metadata_refSeq_file = args.gencode_metadata_refSeq
    MANE_summary_txt = args.MANE_summary_txt
    refseq_gtf_file = args.refseq_gtf_file

    variants = read_file_and_get_variants(UID2_input)
    transcript_file_set = parse_appris_transcript(appris_MANE_transcript_file)
    nm_to_enst =  parse_gencode_nm_to_enst(gencode_metadata_refSeq_file)
    symbol_to_MANEtrx = parse_mane_summary_txt(MANE_summary_txt, 'MANE Select')
    symbol_to_MANEtrx_clinical = parse_mane_summary_txt(MANE_summary_txt, 'MANE Plus Clinical')
    refseq_gtf_dict = parse_gtf_to_transcript_dict(refseq_gtf_file, ['seqnames', 'start', 'end', 'type', 'transcript_id', 'exon_number'])
    get_representative_transcript(variants, json_file_path, transcript_file_set, output_file, nm_to_enst, symbol_to_MANEtrx, symbol_to_MANEtrx_clinical, refseq_gtf_dict)
