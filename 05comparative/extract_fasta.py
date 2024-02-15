import argparse
from Bio import SeqIO

def extract_sequences(ids_file, proteome_file, output_file):
    # Read sequence IDs from ids_file
    with open(ids_file, 'r') as f:
        ids = {line.strip() for line in f}

    # Write sequences to output_file
    with open(output_file, 'w') as out_fasta:
        for record in SeqIO.parse(proteome_file, 'fasta'):
            if record.id in ids:
                SeqIO.write(record, out_fasta, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Extract sequences from a proteome file based on provided IDs.')
    parser.add_argument('--ids_file', help='Path to the file containing sequence IDs', required=True)
    parser.add_argument('--proteome_file', help='Path to the proteome file', required=True)
    parser.add_argument('--output_file', help='Path to the output file to save extracted sequences', required=True)
    args = parser.parse_args()

    extract_sequences(args.ids_file, args.proteome_file, args.output_file)
    print("Sequences extracted successfully!")

if __name__ == "__main__":
    main()

