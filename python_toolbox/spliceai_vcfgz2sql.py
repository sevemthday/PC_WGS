import argparse
import csv
import gzip
import sqlite3

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_spliceai_vcfgz',
        type = str, metavar = 'FILE', required = True,
        help = 'the vcf.gz file from spliceai vcf'
    )
    parser.add_argument(
        '--output_spliceai_concise',
        type = str, metavar = 'FILE', required = True,
        help = 'the UID2 and parsed spliceai score in the text format'
    )
    parser.add_argument(
        '--output_spliceai_sqldb',
        type = str, metavar = 'FILE', required = True,
        help = 'the UID2 and parsed spliceai score in the text sql db'
    )
    args = parser.parse_args()
    return args

def get_column_name():
    return ['UID2', 'Alternate_allele', 'Spliceai_transcript',
            'Delta_score_acceptor_gain', 'Delta_score_acceptor_loss', 'Delta_score_donor_gain', 'Delta_score_donor_loss',
            'Delta_position_acceptor_gain', 'Delta_position_acceptor_loss', 'Delta_position_donor_gain', 'Delta_position_donor_loss',
            'high_precision', 'recommended', 'high_recall']

#def get_DS_cutoff_boolean(deltaScore: str, cutoff: float) -> str:
#    """
#    Checks if any score in the deltaScore list is greater than or equal to the cutoff.
#    Returns 'True' as a string if any score meets the condition, otherwise 'False'.
#    """
#    DS_boolean = str(any(list(map(lambda x:float(x)>=cutoff, deltaScore))))
#    return DS_boolean

def get_DS_cutoff_boolean(deltaScore: str, cutoff: float) -> str:
    """
    Checks if any score in the deltaScore list is greater than or equal to the cutoff.
    Returns 'True' as a string if any score meets the condition, otherwise 'False'.
    Handles cases where the score cannot be converted to a float.
    """
    def is_valid_float(x):
        try:
            return float(x) >= cutoff
        except ValueError:
            return False

    DS_boolean = str(any(list(map(is_valid_float, deltaScore))))
    return DS_boolean

def parse_spliceai_score_raw(spliceai_score_raw: str) -> list:
    """
    Parses the raw SpliceAI score string and returns a list of formatted SpliceAI scores.
    
    If the input is '.', returns a list of thirteen '.' strings.
    Otherwise, processes the SpliceAI scores by splitting, replacing '0.00' with '0', 
    and adding boolean flags based on delta score cutoffs.
    """
    if "SpliceAI=" in spliceai_score_raw:
        # Split the raw SpliceAI scores
        spliceai_score_all = spliceai_score_raw.split("SpliceAI=")[1].split(",")

        spliceai_score_merge = []
        for spliceai_score in spliceai_score_all:
            # Split individual SpliceAI scores and replace '0.00' with '0'
            spliceai_score_split = spliceai_score.split("|")
            spliceai_score_split = list(map(lambda x: x.replace('0.00', '0'), spliceai_score_split))

            # Get DS boolean flags for the specified cutoffs
            DS_boolean = list(map(get_DS_cutoff_boolean, [spliceai_score_split[2:6]]*3, [0.8, 0.5, 0.2]))

            # Append the processed scores and boolean flags to the merge list
            spliceai_score_merge.append(spliceai_score_split + DS_boolean)
        
        # Combine and format the processed scores for output
        spliceai_score_merge = [';'.join(x) for x in zip(*spliceai_score_merge)]

        return spliceai_score_merge
    else:
        return ['.'] * 13

def vcfgz2concise(compress_vcf, concise_txt):
    # Open the input VCF file and the output file
    with gzip.open(compress_vcf, 'rb') as f_i, open(concise_txt, 'w') as f_o:
        # Write column names to the output file
        column_name = get_column_name()
        f_o.write('\t'.join(column_name) + '\n')
        
        for row in f_i:
            # Remove leading/trailing whitespace and decode from bytes to string
            row = row.strip().decode('utf-8').split('\t')
            # Process only non-header lines
            if not row[0].startswith('#'):
                if row[4] != '*':
                    # Create a unique identifier for each variant
                    UID2 = f'{row[0]}-{row[1]}-{row[3]}-{row[4]}'
                    
                    # Parse the SpliceAI score from the raw score string in column 8
                    spliceai_score_raw = row[7]
                    spliceai_score_parse = parse_spliceai_score_raw(spliceai_score_raw)
                    
                    # Write the UID and parsed scores to the output file
                    f_o.write('\t'.join([UID2] + spliceai_score_parse) + '\n')



def create_table(cursor, table_name, headers):
    columns = ', '.join(f'{header} TEXT' for header in headers)
    create_table_sql = f'CREATE TABLE IF NOT EXISTS {table_name} ({columns}, PRIMARY KEY (UID2))'
    cursor.execute(create_table_sql)

def insert_data(cursor, table_name, headers, data):
    placeholders = ', '.join(['?' for _ in headers])
    insert_sql = f'INSERT OR IGNORE INTO {table_name} ({", ".join(headers)}) VALUES ({placeholders})'
    cursor.execute(insert_sql, data)

def concise2sqldb(concise_txt, sql_db):
    table_name = 'spliceai_score'

    conn = sqlite3.connect(sql_db)
    cursor = conn.cursor()

    with open(concise_txt, 'r', newline='') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        headers = next(csv_reader)  # Get headers from CSV
        create_table(cursor, table_name, headers)

        for row in csv_reader:
            insert_data(cursor, table_name, headers, row)

    conn.commit()
    conn.close()

def main():
    args = parse_args()
    vcfgz2concise(args.input_spliceai_vcfgz, args.output_spliceai_concise)
    concise2sqldb(args.output_spliceai_concise, args.output_spliceai_sqldb)

if __name__ == "__main__":
    main()
