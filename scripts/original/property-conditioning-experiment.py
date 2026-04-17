#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np
from project.scripts.inference.generate_samples import main as generate_samples
from project.config import load_model_for_inference
from project.data import load_fasta_to_df
from project.sequence_properties import calculate_charge, calculate_hydrophobicity, calculate_length
from hydra import compose, initialize
import torch

def main():
    print("Starting property conditioning experiment...")
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate sequences using property conditioning')
    parser.add_argument('--num_samples', type=int, default=200, help='Number of samples to generate for each property value')
    parser.add_argument('--num_quantiles', type=int, default=10, help='Number of quantiles to use for conditioning')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size for generation')
    parser.add_argument('--amp_sequences_file', type=str, default='data/activity-data/curated-AMPs.fasta', help='Path to the file containing AMP sequences')
    parser.add_argument('--checkpoint_path', type=str, default='models/generative_model.ckpt', help='Path to the model checkpoint')
    parser.add_argument('--output_dir', type=str, default='experiments/results/files/property-conditioning', help='Directory to save output files')
    args = parser.parse_args()
    
    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory created at: {output_dir}")
    
    # Load model
    print("Loading model...")
    with initialize(version_base=None, config_path="../../../config"):
        config = compose(config_name="train")
        model = load_model_for_inference(config, args.checkpoint_path)
    
    # Load AMP sequences
    amp_df = load_fasta_to_df(args.amp_sequences_file)
    amp_sequences = amp_df['Sequence'].tolist()

    # Calculate property values
    print("Calculating property values for all sequences...")
    charge = calculate_charge(amp_sequences)
    hydrophobicity = calculate_hydrophobicity(amp_sequences)
    length = calculate_length(amp_sequences)
    print(f"Property ranges - Charge: [{min(charge):.2f}, {max(charge):.2f}], "
          f"Hydrophobicity: [{min(hydrophobicity):.2f}, {max(hydrophobicity):.2f}], "
          f"Length: [{min(length)}, {max(length)}]")

    # Calculate quantiles
    quantiles = np.round(np.linspace(0, 1, args.num_quantiles) * (len(amp_sequences) - 1)).astype(int)

    print(quantiles)

    # Sort the sequences by the property values
    sorted_charge = np.sort(charge)
    sorted_hydrophobicity = np.sort(hydrophobicity)
    sorted_length = np.sort(length)

    
    #Create dataframes for each property
    charge_df = pd.DataFrame({'Sequence': [], 'Conditioned-Charge': [], 'Obtained-Charge': []})
    hydrophobicity_df = pd.DataFrame({'Sequence': [], 'Conditioned-Hydrophobicity': [], 'Obtained-Hydrophobicity': []})
    length_df = pd.DataFrame({'Sequence': [], 'Conditioned-Length': [], 'Obtained-Length': []})

    # Generate samples for each property value
    print("\nGenerating samples for each property quantile...")
    for i in range(args.num_quantiles):
        print(f"\nProcessing quantile {i+1}/{args.num_quantiles}:")
        
        # Get the sequences for the current quantile
        charge_quantile = sorted_charge[quantiles[i]]
        hydrophobicity_quantile = sorted_hydrophobicity[quantiles[i]]
        length_quantile = sorted_length[quantiles[i]]
        
        print(f"  Generating samples with charge = {charge_quantile:.2f}...")
        charge_samples, _ = generate_samples(mode='PartialConditional', charge=charge_quantile, num_samples=args.num_samples, batch_size=args.batch_size, model=model,
                                             conditioning_output_path=None,
                                             output_fasta=None)
        print(f"  Generated {len(charge_samples)} samples for charge = {charge_quantile:.2f}")
        
        print(f"  Generating samples with hydrophobicity = {hydrophobicity_quantile:.2f}...")
        hydrophobicity_samples, _ = generate_samples(mode='PartialConditional', hydrophobicity=hydrophobicity_quantile, num_samples=args.num_samples, batch_size=args.batch_size, model=model,
                                                     conditioning_output_path=None,
                                                     output_fasta=None)
        print(f"  Generated {len(hydrophobicity_samples)} samples for hydrophobicity = {hydrophobicity_quantile:.2f}")
        
        print(f"  Generating samples with length = {length_quantile}...")
        length_samples, _ = generate_samples(mode='PartialConditional', length=length_quantile, num_samples=args.num_samples, batch_size=args.batch_size, model=model, 
                                             conditioning_output_path=None,
                                             output_fasta=None)
        print(f"  Generated {len(length_samples)} samples for length = {length_quantile}")

        # Compute properties of the generated samples and update the dataframes
        print("  Calculating properties of generated samples and updating dataframes...")
        charge_df = pd.concat([charge_df, pd.DataFrame({'Sequence': charge_samples, 'Conditioned-Charge': [charge_quantile] * len(charge_samples), 'Obtained-Charge': calculate_charge(charge_samples)})])
        hydrophobicity_df = pd.concat([hydrophobicity_df, pd.DataFrame({'Sequence': hydrophobicity_samples, 'Conditioned-Hydrophobicity': [hydrophobicity_quantile] * len(hydrophobicity_samples), 'Obtained-Hydrophobicity': calculate_hydrophobicity(hydrophobicity_samples)})])
        length_df = pd.concat([length_df, pd.DataFrame({'Sequence': length_samples, 'Conditioned-Length': [length_quantile] * len(length_samples), 'Obtained-Length': calculate_length(length_samples)})])

    # Save the dataframes to csv files
    print("\nSaving results to files...")
    charge_df.to_csv(os.path.join(output_dir, 'charge_df.csv'), index=False)
    print(f"Charge conditioned samples saved to {os.path.join(output_dir, 'charge_df.csv')}")
    
    hydrophobicity_df.to_csv(os.path.join(output_dir, 'hydrophobicity_df.csv'), index=False)
    print(f"Hydrophobicity conditioned samples saved to {os.path.join(output_dir, 'hydrophobicity_df.csv')}")
    
    length_df.to_csv(os.path.join(output_dir, 'length_df.csv'), index=False)
    print(f"Length conditioned samples saved to {os.path.join(output_dir, 'length_df.csv')}")
    
    print("\nExperiment completed successfully!")

if __name__ == "__main__":
    main()
