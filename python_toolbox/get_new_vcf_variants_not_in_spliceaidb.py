import argparse
import gzip
import sqlite3

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_vcfgz',
        type = str, metavar = 'FILE', required = True,
        help = 'the vcf.gz file from spliceai vcf'
    )
    parser.add_argument(
        '--output_vcf',
        type = str, metavar = 'FILE', required = True,
        help = 'the filtered vcf'
    )
    parser.add_argument(
        '--spliceai_sqldb',
        type = str, metavar = 'FILE', required = True,
        help = 'the spliceai sql database'
    )
    args = parser.parse_args()
    return args

def get_existing_UID2s(cursor, table_name):
    cursor.execute(f'SELECT UID2 FROM {table_name}')
    existing_uids = set(row[0] for row in cursor.fetchall())
    return existing_uids

def main():
    args = parse_args()

    input_vcf = args.input_vcfgz
    db_name = args.spliceai_sqldb
    table_name = 'spliceai_score'
    output_vcf = args.output_vcf

    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    existing_uids = get_existing_UID2s(cursor, table_name)
    conn.close()

    with gzip.open(input_vcf, 'rb') as f_i, open(output_vcf, 'w') as f_o:
        for row in f_i:
            row = row.decode('utf-8')
            if row[0] == '#':
                f_o.write(row)
            else:
                row_split = row.strip().split('\t')
                UID2 = f'{row_split[0]}-{row_split[1]}-{row_split[3]}-{row_split[4]}'
                if UID2 not in existing_uids:
                    if row_split[4] != '*':
                        f_o.write(row)

if __name__ == "__main__":
    main()
