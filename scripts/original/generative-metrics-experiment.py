#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
import random
from Bio import SeqIO

def load_and_sample_fasta(fasta_path, max_samples=1000):
    """Load sequences from a FASTA file and sample a subset if needed."""
    sequences = []
    with open(fasta_path, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences.append(record)
    
    if len(sequences) > max_samples:
        sampled_sequences = random.sample(sequences, max_samples)
        return sampled_sequences
    return sequences

def write_sampled_fasta(sequences, output_path):
    """Write a list of SeqIO records to a FASTA file."""
    with open(output_path, 'w') as output_file:
        SeqIO.write(sequences, output_file, "fasta")

def run_metrics_for_model(model_name, fasta_path, output_dir, no_generated_samples=1000):
    """Run metrics evaluation for a specific model and save results to CSV."""
    
    # Create output filename without timestamp
    output_csv = os.path.join(output_dir, f"{model_name}-metrics.csv")
    
    # Construct command with no_generated_samples parameter
    cmd = [
        "python", "project/scripts/evaluation/get_generative_model_metrics.py",
        "--path_to_fasta", fasta_path,
        "--output_csv", output_csv,
        "--no_generated_samples", str(no_generated_samples)
    ]
    
    print(f"Running metrics for {model_name} with {no_generated_samples} samples...")
    subprocess.run(cmd, check=True)
    print(f"Metrics for {model_name} saved to {output_csv}")
    
    return output_csv

def main():
    # Create output directory for results
    output_dir = "experiments/results/files/generative-metrics"
    os.makedirs(output_dir, exist_ok=True)
    
    # Number of samples to use for evaluation
    no_generated_samples = 50000
    
    # Define models and their data sources
    models_and_generated_data = [
        ("OmegAMP-PropertyConditional", "experiments/results/files/generative-samples/OmegAMP/property-conditional-samples.fasta"),
        ("OmegAMP-AllScales-pI", "experiments/results/files/generative-samples/OmegAMP-ablations/all-scales-but-pI-generated-samples.fasta"),
        ("OmegAMP-AllScales-WW", "experiments/results/files/generative-samples/OmegAMP-ablations/all-scales-but-WW-generated-samples.fasta"),
        ("HQ-Non-AMPs", "data/activity-data/curated-Non-AMPs.fasta"),
        ("HQ-AMPs", "data/activity-data/curated-AMPs.fasta"),
        ("Diff-AMP", "experiments/results/files/generative-samples/Diff-AMP/diff-amp-samples.fasta"),
        ("OmegAMP-OneHot", "experiments/results/files/generative-samples/OmegAMP-ablations/one-hot-generated-samples.fasta"),
        ("OmegAMP-Numeric", "experiments/results/files/generative-samples/OmegAMP-ablations/numeric-generated-samples.fasta"),
        ("OmegAMP-Unconditional", "experiments/results/files/generative-samples/OmegAMP/omegamp-generated-samples.fasta"),
        ("OmegAMP-SubsetConditional", "experiments/results/files/generative-samples/OmegAMP/subset-hq-conditional-samples.fasta"),
        ("HydrAMP", "experiments/results/files/generative-samples/HydrAMP/hydramp-generated-samples.fasta"),
        ("AMPGAN", "experiments/results/files/generative-samples/amp-gan/amp-gan-samples.fasta"),
        ("AMP-Diffusion", "experiments/results/files/generative-samples/AMP-Diffusion/amp-diffusion.fasta"),
    ]
    
    # Run metrics for each model and collect CSV paths
    csv_files = []
    for model_name, fasta_path in models_and_generated_data:
        if os.path.exists(fasta_path):
            csv_path = run_metrics_for_model(model_name, fasta_path, output_dir, no_generated_samples)
            csv_files.append((model_name, csv_path))
        else:
            print(f"Warning: File not found for {model_name}: {fasta_path}")
    # Combine all results into a single comparison CSV
    if csv_files:
        combined_df = pd.DataFrame()
        for model_name, csv_path in csv_files:
            model_df = pd.read_csv(csv_path)
            model_df['Model'] = model_name
            combined_df = pd.concat([combined_df, model_df], ignore_index=True)
        
        # Save combined results
        combined_csv = os.path.join(output_dir, f"all-models-comparison.csv")
        combined_df.to_csv(combined_csv, index=False)
        print(f"Combined metrics for all models saved to {combined_csv}")

if __name__ == "__main__":
    main()
