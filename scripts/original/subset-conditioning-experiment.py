import os
import argparse
from project.scripts.inference.generate_samples import main as generate_samples
from project.scripts.inference.predict_sequences import main as predict_sequences
from project.config import load_model_for_inference
from project.constants import AMPs_FILE, HQ_AMPs_FILE
from hydra import compose, initialize


species_name = {
    "acinetobacterbaumannii": "A. Baumannii",
    "escherichiacoli": "E. Coli",
    "klebsiellapneumoniae": "K. Pneumoniae",
    "pseudomonasaeruginosa": "P. Aeruginosa",
    "staphylococcusaureus": "S. Aureus"
}

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate sequences using subset conditioning with different subsets and predict their properties.')
    parser.add_argument('--num_samples', type=int, default=2000, help='Number of samples to generate for each conditioning')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size for generation')
    parser.add_argument('--checkpoint_path', type=str, default='models/generative_model.ckpt', help='Path to the model checkpoint')
    parser.add_argument('--output_dir', type=str, default='experiments/results/files/subset-conditioning', help='Directory to save output files')
    parser.add_argument('--skip_generation', action='store_true', help='Skip sequence generation and only run predictions (assumes files exist)')
    args = parser.parse_args()
    
    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    with initialize(version_base=None, config_path="../../../config"):
        config = compose(config_name="train")
        model = load_model_for_inference(config, args.checkpoint_path)
    
    # Define species files (assuming they follow a standard naming convention)
    species_files = {
        species: f"data/activity-data/strain-species-data/species/{species}_positive.fasta" 
        for species in species_name.keys()
    }
    
    # --- All AMPs Conditioning ---
    amp_output = f"{output_dir}/all-amps-samples.fasta"
    amp_conditioning = f"{output_dir}/all-amps-conditioning.pt"

    if not args.skip_generation:
        # Generate samples using regular AMPs file
        print(f"Generating {args.num_samples} samples conditioned on all AMPs...")
        generate_samples(
            mode='Unconditional',
            checkpoint_path=args.checkpoint_path,
            output_fasta=amp_output,
            num_samples=args.num_samples,
            batch_size=args.batch_size,
            model=model
        )

    # Run prediction on generated samples
    print(f"Running predictions on all AMPs samples...")
    predict_sequences(
        fasta_file=amp_output,
        classifier_choice="all",
        output_csv=f"{output_dir}/all-amps-predictions.csv",
        predict_proba=True
    )
    
    # --- High-Quality AMPs Conditioning ---
    hq_amp_output = f"{output_dir}/hq-amps-samples.fasta"
    hq_amp_conditioning = f"{output_dir}/hq-amps-conditioning.pt"

    if not args.skip_generation:
        # Generate samples using high-quality AMPs file
        print(f"Generating {args.num_samples} samples conditioned on high-quality AMPs...")
        generate_samples(
            mode='SubsetConditional',
            subset_sequences=HQ_AMPs_FILE,
            checkpoint_path=args.checkpoint_path,
            output_fasta=hq_amp_output,
            conditioning_output_path=hq_amp_conditioning,
            num_samples=args.num_samples,
            batch_size=args.batch_size,
            model=model
        )

    # Run prediction on generated samples
    print(f"Running predictions on high-quality AMPs samples...")
    predict_sequences(
        fasta_file=hq_amp_output,
        classifier_choice="all",
        output_csv=f"{output_dir}/hq-amps-predictions.csv",
        predict_proba=True
    )
    
    # --- Species-Specific Conditioning ---
    for species_id, species_display in species_name.items():
        species_output = f"{output_dir}/{species_id}-samples.fasta"
        species_conditioning = f"{output_dir}/{species_id}-conditioning.pt"

        if not args.skip_generation:
            print(f"Generating {args.num_samples} samples conditioned on {species_display} AMPs...")
            generate_samples(
                mode='SubsetConditional',
                subset_sequences=species_files[species_id],
                checkpoint_path=args.checkpoint_path,
                output_fasta=species_output,
                conditioning_output_path=species_conditioning,
                num_samples=args.num_samples,
                batch_size=args.batch_size,
                model=model
            )

        # Run prediction on generated samples
        print(f"Running predictions on {species_display} samples...")
        predict_sequences(
            fasta_file=species_output,
            classifier_choice="all",
            output_csv=f"{output_dir}/{species_id}-predictions.csv",
            predict_proba=True
        )

if __name__ == "__main__":
    main()
