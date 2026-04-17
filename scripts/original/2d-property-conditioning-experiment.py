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
from project.classifiers import AMPClassifier
from project.constants import CLASSIFIER_MODELS

def main():
    print("Starting 2D property conditioning experiment...")
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate sequences conditioning on charge and hydrophobicity, optionally predict AMP probability on existing file.')
    parser.add_argument('--num_samples', type=int, default=2048, help='Number of samples to generate for each property value combination')
    parser.add_argument('--batch_size', type=int, default=512, help='Batch size for generation')
    parser.add_argument('--min_length', type=int, default=10, help='Minimum length for generated sequences')
    parser.add_argument('--max_length', type=int, default=30, help='Maximum length for generated sequences')
    parser.add_argument('--amp_sequences_file', type=str, default='data/activity-data/curated-AMPs.fasta', help='Path to the file containing AMP sequences')
    parser.add_argument('--checkpoint_path', type=str, default='models/generative_model.ckpt', help='Path to the model checkpoint')
    parser.add_argument('--output_dir', type=str, default='experiments/results/files/2d-property-conditioning', help='Directory to save output files')
    parser.add_argument('--predict_only', action='store_true', help='If set, skip generation and run only AMP prediction on the default output CSV file.')
    args = parser.parse_args()
    
    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory created at: {output_dir}")
    
    # Define the standard output file path (used in both modes)
    output_filename = 'charge_hydrophobicity_conditioned_samples.csv'
    output_path = os.path.join(output_dir, output_filename)

    # --- Check for Prediction-Only Mode ---
    if args.predict_only:
        print(f"\n--- Running in Prediction-Only Mode ---")
        print(f"Attempting to load existing data from: {output_path}")
        if not os.path.exists(output_path):
            print(f"Error: File not found at {output_path}")
            print("Please run the script without --predict_only first to generate the file.")
            return # Exit if file doesn't exist
        results_df = pd.read_csv(output_path)
        print(f"Loaded {len(results_df)} sequences.")
        # output_path is already defined and will be used for saving later
    else:
        print("\n--- Running Full Generation and Prediction Mode ---")
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
        # Length calculation is still useful for context, even if not directly used for conditioning quantiles
        length = calculate_length(amp_sequences) 
        print(f"Property ranges - Charge: [{min(charge):.2f}, {max(charge):.2f}], "
              f"Hydrophobicity: [{min(hydrophobicity):.2f}, {max(hydrophobicity):.2f}], "
              f"Length: [{min(length)}, {max(length)}]")

        # Calculate quantiles (Min, Median, Max)
        num_quantiles = 3 # Fixed to 3 for Q1, Q2, Q3 (Min, Median, Max)
        quantile_indices = np.round(np.linspace(0.25, 0.75, num_quantiles) * (len(amp_sequences) - 1)).astype(int)

        print(f"Using quantile indices: {quantile_indices}")

        # Sort the sequences by the property values
        sorted_charge = np.sort(charge)
        sorted_hydrophobicity = np.sort(hydrophobicity)
        # sorted_length is not needed for conditioning grid

        # Get the specific quantile values for charge and hydrophobicity
        charge_quantile_values = sorted_charge[quantile_indices]
        hydrophobicity_quantile_values = sorted_hydrophobicity[quantile_indices]

        print(f"Charge Quantile Values (Q1, Q2, Q3): {np.round(charge_quantile_values, 2)}")
        print(f"Hydrophobicity Quantile Values (Q1, Q2, Q3): {np.round(hydrophobicity_quantile_values, 2)}")
        
        # Create a list to store results
        results_list = []

        # Generate samples for each property value combination (3x3 grid)
        print(f"\nGenerating samples for {num_quantiles}x{num_quantiles} property quantile grid...")
        for i, charge_val in enumerate(charge_quantile_values):
            for j, hydro_val in enumerate(hydrophobicity_quantile_values):
                print(f"\nProcessing combination Charge Q{i+1} ({charge_val:.2f}), Hydrophobicity Q{j+1} ({hydro_val:.2f}):")
                
                print(f"  Generating samples with charge = {charge_val:.2f} and hydrophobicity = {hydro_val:.2f}...")
                # Generate samples conditioning on charge and hydrophobicity, length is unconditioned
                generated_samples, _ = generate_samples(
                    mode='PartialConditional', 
                    charge=f"{charge_val}", 
                    hydrophobicity=f"{hydro_val}", 
                    length=f"{args.min_length}:{args.max_length}", # Length range for generation, not filtering
                    num_samples=args.num_samples, 
                    batch_size=args.batch_size, 
                    model=model,
                    conditioning_output_path=None,
                    output_fasta=None
                )
                print(f"  Generated {len(generated_samples)} raw samples.")

                # Compute properties of the generated samples
                print("  Calculating properties of generated samples...")
                obtained_charge = calculate_charge(generated_samples)
                obtained_hydrophobicity = calculate_hydrophobicity(generated_samples)
                obtained_length = calculate_length(generated_samples)

                # Store results for all generated samples (no length filtering)
                print(f"  Storing results for all {len(generated_samples)} generated samples...")
                for k in range(len(generated_samples)):
                    results_list.append({
                        'Sequence': generated_samples[k],
                        'Conditioned-Charge': charge_val,
                        'Conditioned-Hydrophobicity': hydro_val,
                        'Obtained-Charge': obtained_charge[k],
                        'Obtained-Hydrophobicity': obtained_hydrophobicity[k],
                        'Obtained-Length': obtained_length[k]
                    })
                print(f"  Stored {len(generated_samples)} samples.")


        # Convert results list to DataFrame
        print("\nConverting results to DataFrame...")
        results_df = pd.DataFrame(results_list)
        # output_path is already defined
    # --- End of Mode Check ---

    # --- Predict AMP probabilities (Runs in both modes) ---
    print("\nPredicting AMP probabilities...")
    # Load the broad AMP classifier
    classifier_name = 'broad-classifier'
    if classifier_name not in CLASSIFIER_MODELS:
        raise ValueError(f"Classifier '{classifier_name}' not found in CLASSIFIER_MODELS.")
    
    model_path = CLASSIFIER_MODELS[classifier_name]
    print(f"Loading classifier: {classifier_name} from {model_path}")
    amp_classifier = AMPClassifier(model_path=model_path)
    amp_classifier.eval() # Set model to evaluation mode

    # Get sequences from the DataFrame
    sequences_to_predict = results_df['Sequence'].tolist()
    
    # Predict probabilities
    print(f"Running prediction for {len(sequences_to_predict)} sequences...")
    amp_probabilities = amp_classifier.predict_proba(sequences_to_predict)
    
    # Add probabilities to the DataFrame
    results_df['amp-prediction'] = amp_probabilities
    print("AMP probabilities added to the DataFrame.")
    # --- End of prediction step ---

    # Save the dataframe to a csv file
    # output_filename and output_path are defined earlier
    print(f"\nSaving results to {output_path}...")
    results_df.to_csv(output_path, index=False)
    # Modify the final message based on the mode
    if args.predict_only:
        print(f"Updated file '{output_path}' with AMP predictions saved.")
    else:
        print(f"Generated samples and predictions saved to '{output_path}'.")
    
    print("\nExperiment completed successfully!")

if __name__ == "__main__":
    main()
