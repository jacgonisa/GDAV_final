#1. I wrote extract_fasta.py

______
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
________


python 05comparative/extract_fasta.py --ids_file 04de_analysis/results/degs_ids.txt --proteome_file proteome.faa --output_file 05comparative/degs.faa


#2. I wrote divide_fasta.py
_____________________________________
from Bio import SeqIO
import os
import argparse

def divide_fasta(input_file, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over sequences in the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        gene_name = record.id.split()[0]  # Extract the gene name from the header line
        output_file = os.path.join(output_dir, f"{gene_name}.fa")   # Output file path with gene name

        # Write the sequence to the output file
        SeqIO.write(record, output_file, "fasta")

        print(f"Sequence '{record.id}' saved to '{output_file}'")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Divide a FASTA file into individual files based on gene names.")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_dir", help="Output directory for individual FASTA files")

    # Parse command-line arguments
    args = parser.parse_args()

    # Divide the FASTA file into individual files
    divide_fasta(args.input_file, args.output_dir)

_______________________________________

python 05comparative/divide_fasta.py 05comparative/degs.faa 05comparative/

Sequence 'AQUIFEX_01423' saved to '05comparative/AQUIFEX_01423.fa'
Sequence 'AQUIFEX_01723' saved to '05comparative/AQUIFEX_01723.fa'
Sequence 'AQUIFEX_01749' saved to '05comparative/AQUIFEX_01749.fa'
Sequence 'AQUIFEX_01759' saved to '05comparative/AQUIFEX_01759.fa'
Sequence 'AQUIFEX_01761' saved to '05comparative/AQUIFEX_01761.fa'


