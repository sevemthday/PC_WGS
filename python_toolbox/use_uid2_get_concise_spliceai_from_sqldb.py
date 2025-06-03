import argparse
import sqlite3

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_UID2',
        type = str, metavar = 'FILE', required = True,
        help = 'the UID2 list'
    )
    parser.add_argument(
        '--output_concise_spliceai',
        type = str, metavar = 'FILE', required = True,
        help = 'the output txt of concise spliceai score'
    )
    parser.add_argument(
        '--spliceai_sqldb',
        type = str, metavar = 'FILE', required = True,
        help = 'the path of spliceai sql database'
    )
    args = parser.parse_args()
    return args

def get_column_name():
    return ['UID2', 'Alternate_allele', 'Spliceai_transcript',
            'Delta_score_acceptor_gain', 'Delta_score_acceptor_loss', 'Delta_score_donor_gain', 'Delta_score_donor_loss',
            'Delta_position_acceptor_gain', 'Delta_position_acceptor_loss', 'Delta_position_donor_gain', 'Delta_position_donor_loss',
            'high_precision', 'recommended', 'high_recall']

def read_UID2(input_UID2: str) -> list:
    UID2s = []
    with open(input_UID2, 'r') as f:
        for row in f:
            row =  row.strip()
            UID2s.append(row)

    return UID2s

def write_rows(rows: list, output_file: str):
    column_name = get_column_name()

    with open(output_file, 'w') as f:

        line = '\t'.join(column_name)
        f.write(line + '\n')

        for row in rows:
            line = '\t'.join(row)
            f.write(line + '\n')

# def get_concise_spliceai_table(UID2: list, spliceai_sqldb: str) -> list:

#     conn = sqlite3.connect(spliceai_sqldb)
#     cursor = conn.cursor()

#     table_name = 'spliceai_score'

#     query = f"SELECT * FROM {table_name} WHERE UID2 IN ({','.join(['?' for _ in UID2])})"
#     cursor.execute(query, UID2)
#     rows = cursor.fetchall()

#     cursor.close()
#     conn.close()

#     return rows

def get_concise_spliceai_table(UID2: list, spliceai_sqldb: str, chunk_size: int = 500) -> list:
    conn = sqlite3.connect(spliceai_sqldb)
    cursor = conn.cursor()

    table_name = 'spliceai_score'
    all_rows = []

    # Split the UID2 list into chunks
    for i in range(0, len(UID2), chunk_size):
        chunk = UID2[i:i + chunk_size]
        query = f"SELECT * FROM {table_name} WHERE UID2 IN ({','.join(['?' for _ in chunk])})"
        cursor.execute(query, chunk)
        all_rows.extend(cursor.fetchall())

    cursor.close()
    conn.close()

    return all_rows

def main():
    args = parse_args()
    
    UID2s = read_UID2(args.input_UID2)
    concise_spliceai_table = get_concise_spliceai_table(UID2s, args.spliceai_sqldb)
    write_rows(concise_spliceai_table, args.output_concise_spliceai)

if __name__ == "__main__":
    main()
