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


