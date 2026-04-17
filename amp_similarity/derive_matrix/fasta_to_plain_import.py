from Bio import SeqIO
from io import StringIO

def fasta_to_plain_file(input_path, output_path):
    """Convert FASTA file to plain text file (sequences only, no headers)."""
    with open(output_path, "w") as out_file:
        for record in SeqIO.parse(input_path, "fasta"):
            out_file.write(str(record.seq) + "\n")
    print(f"[INFO] Written sequences (one per line) to: {output_path}")


def fasta_str_to_plain_str(fasta_str):
    """
    Convert FASTA-formatted string to plain string (one sequence per line, no headers).
    """
    records = SeqIO.parse(StringIO(fasta_str), "fasta")
    sequences = [str(rec.seq) for rec in records]
    return "\n".join(sequences) + "\n"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert FASTA to plain format (no headers).")
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file path")
    parser.add_argument("--output", "-o", required=True, help="Output file path")

    args = parser.parse_args()
    fasta_to_plain_file(args.input, args.output)
