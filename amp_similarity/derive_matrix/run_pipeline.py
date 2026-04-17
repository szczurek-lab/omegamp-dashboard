import os
import argparse
import pickle
from collections import Counter, defaultdict
from io import StringIO
from Bio.SeqIO import FastaIO
from Bio import SeqIO

from extract_counts import run_blosum_from_seqs
from run_mafft import run_mafft
from extract_blocks_import import extract_blocks_from_fasta_string
from fasta_to_plain_import import fasta_str_to_plain_str 

def parse_clstr(clstr_file):
    clusters = defaultdict(list)
    with open(clstr_file) as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_cluster = int(line.strip().split()[1])
            else:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                seq_name = parts[2].strip('>').split('...')[0]
                clusters[current_cluster].append(seq_name)
    return clusters

def dict_to_fasta_string(seq_dict):
    """Convert dict of SeqRecords to fasta formatted string"""
    output = StringIO()
    writer = FastaIO.FastaWriter(output, wrap=None)
    writer.write_file(seq_dict.values())
    return output.getvalue()

def merge_counters(counters):
    result = Counter()
    for c in counters:
        result.update(c)
    return result

def load_fasta_to_dict(fasta_file):
    print("[INFO] Loading entire FASTA into memory...")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return seq_dict

def main(fasta_file, clstr_file, output_dir, final_suffix):
    os.makedirs(output_dir, exist_ok=True)
    clusters = parse_clstr(clstr_file)
    final_counts = Counter()
    all_seqs = load_fasta_to_dict(fasta_file)

    for cluster_num, seq_ids in clusters.items():

        if cluster_num == 2287:
            continue  # Skip cluster 2287 it is handled separately because it is too large

        if len(seq_ids) <= 1:
            print(f"Skipping cluster {cluster_num} with only one sequence.")
            continue

        print(f"[INFO] Processing cluster {cluster_num} with {len(seq_ids)} sequences")

        cluster_seqs = {sid: all_seqs[sid] for sid in seq_ids if sid in all_seqs}
        if not cluster_seqs:
            print(f"[WARNING] No sequences found in FASTA for cluster {cluster_num}")
            continue

        cluster_fasta_str = dict_to_fasta_string(cluster_seqs)

        aligned_fasta_str = run_mafft(cluster_fasta_str)

        blocks = extract_blocks_from_fasta_string(aligned_fasta_str)

        cluster_counts = Counter()

        for i, block_fasta_str in enumerate(blocks, start=1):
            plain_str = fasta_str_to_plain_str(block_fasta_str)
            counts = run_blosum_from_seqs(plain_str, 100)
            cluster_counts.update(counts)

        final_counts.update(cluster_counts)

    print("[INFO] Final aggregated counts: DONE")

    counts_file = os.path.join(output_dir, f"final_counts_{final_suffix}.pkl")
    with open(counts_file, "wb") as f:
        pickle.dump(dict(final_counts), f)
    print(f"[INFO] Final counts saved to {counts_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Master pipeline for cluster analysis")
    parser.add_argument("--fasta", "-f", required=True, help="Input FASTA file")
    parser.add_argument("--clstr", "-c", required=True, help="Input .clstr file")
    parser.add_argument("--output", "-o", default="output", help="Output directory")
    parser.add_argument("--final_sufix", "-n", type=str, default="", help="Suffix for final counter file")
    args = parser.parse_args()
    main(args.fasta, args.clstr, args.output, args.final_sufix)
