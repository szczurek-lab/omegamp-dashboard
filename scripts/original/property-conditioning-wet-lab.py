import os
import argparse
import itertools
from project.scripts.inference.generate_samples import main as generate_samples
from project.scripts.inference.predict_sequences import main as predict_sequences
from project.config import load_model_for_inference
from hydra import compose, initialize


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate sequences using property conditioning')
    parser.add_argument('--num_samples', type=int, default=64, help='Number of samples to generate for each conditioning')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size for generation')
    parser.add_argument('--checkpoint_path', type=str, default='models/generative_model.ckpt', help='Path to the model checkpoint')
    parser.add_argument('--output_dir', type=str, default='experiments/results/files/property-conditioning-neurips', help='Directory to save output files')
    args = parser.parse_args()
    
    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    with initialize(version_base=None, config_path="../../../config"):
        config = compose(config_name="train")
        model = load_model_for_inference(config, args.checkpoint_path)

    # single property conditioning
    property_config = {
        'charge':  2-10,
        'hydrophobicity': -0.5-0.8,
        'length':  10-30,
    }
    output_fasta = f"{output_dir}/propertyconditioning-samples.fasta"
    output_conditioning = f"{output_dir}/propertyconditioning-conditions.pt"
    generate_samples(
        mode='PartialConditional',
        subset_sequences='-',
        checkpoint_path=args.checkpoint_path,
        output_fasta=output_fasta,
        conditioning_output_path=output_conditioning,
        num_samples=args.num_samples,
        batch_size=args.batch_size,
        length=property_config['length'],
        charge=property_config['charge'],
        hydrophobicity=property_config['hydrophobicity'],
        guidance_strength=1,
        analog='',
        template='',
        model=model
        )

    # Run prediction on generated samples
    print(f"Running predictions on all AMPs samples...")
    predict_sequences(
        fasta_file=output_fasta,
        classifier_choice="all",
        output_csv=f"{output_dir}/propertyconditioning-predictions.csv",
        predict_proba=True
    )


if __name__ == "__main__":
    main()
