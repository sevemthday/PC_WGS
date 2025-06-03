import argparse
import csv
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        type = str, metavar = 'FILE', required = True,
        help = 'the input file with spliceai annotatio'
    )
    parser.add_argument(
        '--output',
        type = str, metavar = 'FILE', required = True,
        help = 'the output file'
    )
    args = parser.parse_args()
    return args

def read_data(file_path: str) -> list:
    data = []
    with open(file_path, newline='') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter='\t')
        for row in csv_reader:
            data.append(row)
    return data

def process_transcript_string(transcript_string: str, sep: str) -> list:
    transcript_string_split = transcript_string.split(sep)
    transcript_string_head = [i.split('.')[0] for i in transcript_string_split]
    return transcript_string_head

def get_selected_spliceai_transcript_index(row: dict) -> list:
    # get spliceai transcript
    spliceai_trx = process_transcript_string(row['Spliceai_transcript'], ';')
    # get representive trnascirpt info: mane select or not
    is_rTrs_mane_select = row['Mane_Select']
    # get representive trnascirpt info: mane_select, mane_plus_clinical, refseq-paired enst, matching hgsv enst
    trx_mane_select = [row['t_mane_select_ensembl']]
    trx_mane_plus_clinical = [row['t_mane_plus_clinical_ensembl']]
    trx_refseq_paired_enst = process_transcript_string(row['rTrx_paired_enst'], ',')
    trx_matching_rTrx_enst = process_transcript_string(row['matching_rTrx_enst'], ',')
    all_trx = trx_mane_select + trx_mane_plus_clinical + trx_refseq_paired_enst + trx_matching_rTrx_enst
    all_trx = list(set(all_trx))
    # check
    if is_rTrs_mane_select == 'True':
        indices = [index for index, value in enumerate(spliceai_trx) if value in trx_mane_select]
        if indices:
            return indices
        else:
            indices = [index for index, value in enumerate(spliceai_trx) if value in all_trx]
            if indices:
                return indices
            else:
                return ['unknown']
    elif is_rTrs_mane_select == 'False':
        indices = [index for index, value in enumerate(spliceai_trx) if value in all_trx]
        if indices:
            return indices
        else:
            return ['unknown']
    elif is_rTrs_mane_select == 'NA':
        return ['unknown']

def get_spliceai_score(row: dict, selected_transcript_index: list) -> str:
    #Delta_score_acceptor_gain	Delta_score_acceptor_loss	Delta_score_donor_gain	Delta_score_donor_loss	high_precision	recommended	high_recall
    spliceai_score_list = [
         row["Delta_score_acceptor_gain"],
         row["Delta_score_acceptor_loss"],
         row["Delta_score_donor_gain"],
         row["Delta_score_donor_loss"]
    ]
    spliceai_score_array = np.array([col.split(";") for col in spliceai_score_list], dtype=float)
    spliceai_score_max_array = np.max(spliceai_score_array[:, selected_transcript_index], axis=1)
    spliceai_score_max = np.max(spliceai_score_max_array)
    max_index = np.argmax(spliceai_score_max_array)
    max_splcing_event = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"][max_index]
    if float(spliceai_score_max) >= 0.2:
        spliceai_score_report = f"{max_splcing_event}={spliceai_score_max}"
    else:
        spliceai_score_report = '.'
    return spliceai_score_report

def process_data(input: list, output: str):
    data = []
    for row in input:
        if 'T' in row['high_recall']:
            selected_trx_index = get_selected_spliceai_transcript_index(row)
            if selected_trx_index[0] == 'unknown':
                row['spliceai_score_report'] = '.'
            else:
                row['spliceai_score_report'] = get_spliceai_score(row, selected_trx_index)
        else:
            row['spliceai_score_report'] = '.'
        data.append(row)
    df = pd.DataFrame(data)
    df.to_csv(output, sep='\t', index=False)

def main():
    args = parse_args()
    input = args.input
    output = args.output
    data = read_data(input)
    process_data(data, output)

if __name__ == "__main__":
    main()
