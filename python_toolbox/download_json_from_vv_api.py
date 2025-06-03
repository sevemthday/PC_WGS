import argparse
import json
import os
import requests

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
    args = parser.parse_args()
    return args

# 讀取檔案內容並取得variants
def read_file_and_get_variants(UID2_input):
    with open(UID2_input, "r") as f:
        variants = [line.strip() for line in f if line.strip()]
    return variants

def fetch_and_save_data(json_folder, variant_description):
    # API info
    base_path_1 = "https://rest.variantvalidator.org/VariantFormatter/variantformatter"
    base_path_2 = """?content-type=application%2Fjson"""
    genome_build = "GRCh38"
    select_transcripts = "all"
    checkonly = "False"
    headers = {'accept': 'application/json'}
    #
    with open(f'{json_folder}/download_log.txt', "a") as log_file:
        transcript_model = "all" #refseq|ensembl|all
        query = f'{base_path_1}/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}{base_path_2}'
        response = requests.get(query, headers=headers)
        if response.status_code == 200:
            data = response.json()
            # Process the returned data
            print(data)
            with open(f'{json_folder}/{variant_description}.json', "w") as json_file:
                json.dump(data, json_file, indent=4)
        elif response.status_code == 500:
            transcript_model = "refseq" #refseq|ensembl|all
            query = f'{base_path_1}/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}{base_path_2}'
            response = requests.get(query, headers=headers)
            if response.status_code == 200:
                data = response.json()
                # Process the returned data
                print(data)
                with open(f'{json_folder}/{variant_description}.json', "w") as json_file:
                    json.dump(data, json_file, indent=4)
                log_file.write(f"下載 {variant_description} {transcript_model} 資料成功\n")
            else:
                log_file.write(f"下載 {variant_description} {transcript_model} 資料失敗：{response.status_code}\n")
                print(f"請求失敗，狀態碼: {response.status_code}")
        else:
            log_file.write(f"下載 {variant_description} {transcript_model} 資料失敗：{response.status_code}\n")
            print(f"請求失敗，狀態碼: {response.status_code}")

def download_json_from_vv_api(variants, json_folder):
    # get existing downloaded variants from the json_folder.
    exsting_files = os.listdir(json_folder)
    exsting_variants = {f.replace('.json', '') for f in exsting_files if f.endswith('.json')}
    with open(f'{json_folder}/download_log.txt', 'w') as file:
        pass
    for variant_description in variants:
        if variant_description not in exsting_variants:
            fetch_and_save_data(json_folder, variant_description)

if __name__ == '__main__':
    args = parse_args()
    UID2_input = args.UID2_input
    json_folder = args.json_folder
    variants = read_file_and_get_variants(UID2_input)
    download_json_from_vv_api(variants, json_folder)
