import argparse
import csv

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        type = str, metavar = 'FILE', required = True,
        help = 'the annovar input txt'
    )
    parser.add_argument(
        '--output',
        type = str, metavar = 'FILE', required = True,
        help = 'the annovar output txt with UID2'
    )
    args = parser.parse_args()
    return args

def process_line(row):
    row['UID2'] = f"{row['Otherinfo4']}-{row['Otherinfo5']}-{row['Otherinfo7']}-{row['Otherinfo8']}"
    return row

def process_file(input_file, output_file):
    with open(input_file, 'r', newline='', encoding='utf-8') as f_i, open(output_file, 'w', newline='', encoding='utf-8') as f_o:
        reader = csv.DictReader(f_i, delimiter='\t')
        
        # get first data
        first_row = next(reader)
        first_row_processed = process_line(first_row)
        
        # get new column name
        fieldnames = list(first_row_processed.keys())
        
        writer = csv.DictWriter(f_o, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerow(first_row_processed)
        
        # process and write other data
        for row in reader:
            processed_row = process_line(row)
            writer.writerow(processed_row)

def main():
    args = parse_args()
    process_file(args.input, args.output)

if __name__ == "__main__":
    main()
