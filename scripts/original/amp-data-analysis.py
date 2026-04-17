#!/usr/bin/env python3

import os
import pandas as pd
from Bio import SeqIO
import numpy as np
from project.sequence_properties import calculate_length, calculate_charge, calculate_hydrophobicity
from project.metrics import FitnessScore, PseudoPerplexity
from project.wrappers import EsmWrapper

def load_fasta_sequences(fasta_path):
    """Load sequences from a FASTA file."""
    sequences = []
    with open(fasta_path, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences.append(str(record.seq))
    return sequences

def calculate_properties(sequences):
    """Calculate sequence properties for a list of sequences."""
    # Initialize ESM model for pseudo-perplexity calculation
    esm_model = EsmWrapper(320, 100, "cpu")

    properties = {
        "length": calculate_length(sequences).tolist(),
        "charge": calculate_charge(sequences).tolist(),
        "hydrophobicity": calculate_hydrophobicity(sequences, scale="eisenberg").tolist(),
        "fitness_score": FitnessScore()(sequences, all_results=True)[1].tolist(),
        "pseudo_perplexity": PseudoPerplexity(esm_model)(sequences, all_results=True)[1].tolist()
    }
    return properties

def main():
    # Define input directory and create output directory
    input_dir = "experiments/results/files/amp-data"
    output_dir = "experiments/results/files/amp-data-analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each FASTA file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            fasta_path = os.path.join(input_dir, filename)
            dataset_name = os.path.splitext(filename)[0]
            
            print(f"Processing {dataset_name}...")
            
            # Load sequences
            sequences = load_fasta_sequences(fasta_path)
            
            sequences = sequences[:5000] # FIXME


            if not sequences:
                print(f"No sequences found in {fasta_path}")
                continue

            # Calculate properties
            properties = calculate_properties(sequences)

            # Create DataFrame from properties
            df = pd.DataFrame({
                "sequence": sequences,
                "length": properties["length"],
                "charge": properties["charge"],
                "hydrophobicity": properties["hydrophobicity"],
                "fitness_score": properties["fitness_score"],
                "pseudo_perplexity": properties["pseudo_perplexity"]
            })
            
            # Save to CSV
            output_csv = os.path.join(output_dir, f"{dataset_name}-properties.csv")
            df.to_csv(output_csv, index=False)
            print(f"Properties for {dataset_name} saved to {output_csv}")
            
            # Calculate and save summary statistics
            summary_df = pd.DataFrame({
                "property": ["length", "charge", "hydrophobicity", "fitness_score", "pseudo_perplexity"],
                "mean": [
                    np.mean(properties["length"]),
                    np.mean(properties["charge"]),
                    np.mean(properties["hydrophobicity"]),
                    np.mean(properties["fitness_score"]),
                    np.mean(properties["pseudo_perplexity"])
                ],
                "std": [
                    np.std(properties["length"]),
                    np.std(properties["charge"]),
                    np.std(properties["hydrophobicity"]),
                    np.std(properties["fitness_score"]),
                    np.std(properties["pseudo_perplexity"])
                ],
                "min": [
                    np.min(properties["length"]),
                    np.min(properties["charge"]),
                    np.min(properties["hydrophobicity"]),
                    np.min(properties["fitness_score"]),
                    np.min(properties["pseudo_perplexity"])
                ],
                "max": [
                    np.max(properties["length"]),
                    np.max(properties["charge"]),
                    np.max(properties["hydrophobicity"]),
                    np.max(properties["fitness_score"]),
                    np.max(properties["pseudo_perplexity"])
                ]
            })
            
            summary_csv = os.path.join(output_dir, f"{dataset_name}-summary.csv")
            summary_df.to_csv(summary_csv, index=False)
            print(f"Summary statistics for {dataset_name} saved to {summary_csv}")

if __name__ == "__main__":
    main()