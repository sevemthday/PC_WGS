import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_vcfgz',
        type = str, metavar = 'FILE', required = True,
        help = 'the vcf.gz file from spliceai vcf'
    )
    parser.add_argument(
        '--output_UID2',
        type = str, metavar = 'FILE', required = True,
        help = 'the UID2 and parsed spliceai score in the text format'
    )
    args = parser.parse_args()
    return args

def vcfgz2uid2(compress_vcf, uid2_txt):
    # Open the input VCF file and the output file
    with gzip.open(compress_vcf, 'rb') as f_i, open(uid2_txt, 'w') as f_o:
        for row in f_i:
            # Remove leading/trailing whitespace and decode from bytes to string
            row = row.strip().decode('utf-8').split('\t')
            # Process only non-header lines
            if not row[0].startswith('#'):
                if row[4] != '*':
                    # Create a unique identifier for each variant
                    UID2 = f'{row[0]}-{row[1]}-{row[3]}-{row[4]}'
                    
                    # Write the UID and parsed scores to the output file
                    f_o.write(UID2 + '\n')

def main():
    args = parse_args()
    vcfgz2uid2(args.input_vcfgz, args.output_UID2)

if __name__ == "__main__":
    main()