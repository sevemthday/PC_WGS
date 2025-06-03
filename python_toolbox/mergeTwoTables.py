import argparse
import csv
import os
import sqlite3

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_file_1',
        type = str, metavar = 'FILE', required = True,
        help = 'the major table'
    )
    parser.add_argument(
        '--input_file_2',
        type = str, metavar = 'FILE', required = True,
        help = 'the second table'
    )
    parser.add_argument(
        '--output_file',
        type = str, metavar = 'FILE', required = True,
        help = 'the output table'
    )
    parser.add_argument(
        '--temp_sqldb',
        type = str, metavar = 'FILE', required = True,
        help = 'the temp sql database'
    )
    parser.add_argument(
        '--key_column',
        type = str, metavar = 'STRING', required = True,
        help = 'the key column for merge use'
    )
    args = parser.parse_args()
    return args

def left_join_csv(sql_db, file1, file2, output_file, key_column, batch_size=1000):
    
    conn = sqlite3.connect(sql_db)
    cursor = conn.cursor()

    # # Create a table and import data from a CSV file
    def create_table_from_csv(csv_file, table_name):
        with open(csv_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            headers = next(reader)
            columns = ', '.join(f'"{header}" TEXT' for header in headers)
            cursor.execute(f'CREATE TABLE IF NOT EXISTS {table_name} ({columns})')

            for row in reader:
                placeholders = ', '.join('?' * len(row))
                cursor.execute(f'INSERT INTO {table_name} VALUES ({placeholders})', row)

    # Import two CSV files into a SQLite database
    create_table_from_csv(file1, 'table1')
    create_table_from_csv(file2, 'table2')

    # Get second column names
    cursor.execute("PRAGMA table_info(table2);")
    columns = cursor.fetchall()
    filtered_columns = [col[1] for col in columns if col[1] != key_column]
    table2_columns = ', '.join([f"table2.{col}" for col in filtered_columns])
    '''
    columns = [
    (0, 'key_column', 'INTEGER', 1, None, 1),
    (1, 'col1', 'TEXT', 0, None, 0),
    (2, 'col2', 'REAL', 0, None, 0),
    (3, 'col3', 'INTEGER', 0, None, 0)
    ]
    這裡的每個元組包含了以下信息：
    欄位序號
    欄位名稱
    欄位數據類型
    是否允許NULL
    默認值
    是否為主鍵
    '''

    # Perform a LEFT JOIN operation
    query = f"""
    SELECT table1.*, {table2_columns}
    FROM table1
    LEFT JOIN table2
    ON table1.{key_column} = table2.{key_column}
    """
    cursor.execute(query)

    # Export the results to a new CSV file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Write column names
        column_names = [desc[0] for desc in cursor.description]
        writer.writerow(column_names)
        # Write data rows
        #for row in cursor.fetchall():
        #    writer.writerow(row)

        # Fetch query results in batches and write them to a CSV file.
        while True:
            rows = cursor.fetchmany(batch_size)
            if not rows:
                break
            for row in rows:
                writer.writerow(row)

    conn.close()

def main():
    args = parse_args()
    left_join_csv(args.temp_sqldb, args.input_file_1, args.input_file_2, args.output_file, args.key_column)
    os.remove(args.temp_sqldb)

if __name__ == "__main__":
    main()