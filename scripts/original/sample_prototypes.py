#!/usr/bin/env python
"""
Sample prototypes from full curated datasets for initial parameter sweep testing.
Creates smaller sample files to test the pipeline before running on full datasets.
"""

import argparse
import random
from pathlib import Path
from Bio import SeqIO


def sample_sequences(input_fasta, output_fasta, n_samples, seed=42):
    """
    Randomly sample n sequences from a FASTA file.
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output FASTA file
        n_samples: Number of sequences to sample
        seed: Random seed for reproducibility
    """
    random.seed(seed)
    
    # Read all sequences
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    print(f"Input file: {input_fasta}")
    print(f"  Total sequences: {len(records)}")
    
    # Sample
    if n_samples >= len(records):
        print(f"  Warning: Requested {n_samples} samples but only {len(records)} available")
        print(f"  Using all {len(records)} sequences")
        sampled = records
    else:
        sampled = random.sample(records, n_samples)
        print(f"  Sampled: {len(sampled)} sequences")
    
    # Write to output
    output_path = Path(output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    SeqIO.write(sampled, output_fasta, "fasta")
    print(f"  Written to: {output_fasta}")
    
    return len(sampled)


def main():
    parser = argparse.ArgumentParser(
        description='Sample prototype sequences for parameter sweep testing'
    )
    parser.add_argument('--amp_input', type=str, 
                       default='data/activity-data/curated-AMPs.fasta',
                       help='Input FASTA file with AMP sequences')
    parser.add_argument('--nonamp_input', type=str,
                       default='data/activity-data/curated-Non-AMPs.fasta',
                       help='Input FASTA file with Non-AMP sequences')
    parser.add_argument('--amp_output', type=str,
                       default='data/activity-data/curated-AMPs_samples.fasta',
                       help='Output FASTA file for sampled AMPs')
    parser.add_argument('--nonamp_output', type=str,
                       default='data/activity-data/curated-Non-AMPs_samples.fasta',
                       help='Output FASTA file for sampled Non-AMPs')
    parser.add_argument('--n_samples', type=int, default=50,
                       help='Number of sequences to sample from each file')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("SAMPLING PROTOTYPE SEQUENCES")
    print("=" * 60)
    print(f"Samples per file: {args.n_samples}")
    print(f"Random seed: {args.seed}")
    print()
    
    # Sample AMPs (positive prototypes)
    print("Sampling positive prototypes (AMPs)...")
    n_amps = sample_sequences(
        args.amp_input,
        args.amp_output,
        args.n_samples,
        args.seed
    )
    print()
    
    # Sample Non-AMPs (negative prototypes)
    print("Sampling negative prototypes (Non-AMPs)...")
    n_nonamps = sample_sequences(
        args.nonamp_input,
        args.nonamp_output,
        args.n_samples,
        args.seed
    )
    print()
    
    print("=" * 60)
    print("SAMPLING COMPLETE")
    print("=" * 60)
    print(f"Positive prototypes: {n_amps} sequences → {args.amp_output}")
    print(f"Negative prototypes: {n_nonamps} sequences → {args.nonamp_output}")
    print()


if __name__ == '__main__':
    main()