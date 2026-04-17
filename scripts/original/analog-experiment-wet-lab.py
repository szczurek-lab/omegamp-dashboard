import os
import argparse
import itertools
from project.scripts.inference.generate_samples import main as generate_samples
from project.scripts.inference.predict_sequences import main as predict_sequences
from project.config import load_model_for_inference
from project.constants import AMPs_FILE, HQ_AMPs_FILE
from project.data import load_fasta_to_df
from hydra import compose, initialize


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate sequences using subset conditioning with different subsets')
    parser.add_argument('--num_samples', type=int, default=100, help='Number of samples to generate for each conditioning')
    parser.add_argument('--batch_size', type=int, default=100, help='Batch size for generation')
    parser.add_argument('--checkpoint_path', type=str, default='models/generative_model.ckpt', help='Path to the model checkpoint')
    parser.add_argument('--output_dir', type=str, default='experiments/results/files/analog-conditioning-wet-lab', help='Directory to save output files')
    args = parser.parse_args()
    
    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    with initialize(version_base=None, config_path="../../../config"):
        config = compose(config_name="train")
        model = load_model_for_inference(config, args.checkpoint_path)


    negatives = {
        'Prevotellin-2': 'MLNYLYDRDINRYRAIIKALGLRK',
        'Bacteroidin-2': 'MGNLVAIVGRPNVGKSTLFNRFH',    # active against commensals
        'Lachnospirin-1': 'LKQLNRFKYKIIKIKRIIKL',
        'Enterococcin-1': 'SIKLIKTVKVVEKIIVEFLKKIKKFSVI',
        'Mylodonin-2': 'KRKRGLKLATALSLNNKF',
        'Mammuthusin-3': 'KTLKIIRLLF',
        'Mammutin-1': 'WMTIHALKLSLSFKL',
        'Arctodutin-1':  'LLINSIKRLLLGSIL',
        'Hesperelin-2': 'KKWVIKSVVFFK',
        'Equusin-3': 'KLYQRILWRLISEL',
        'Equusin-1': 'FLKLRWSRFARVLL',
    }

    for name, sequence in negatives.items():
        print(f"Generating {args.num_samples} analogs for {name}...")
        output_fasta = f"{output_dir}/{name}-samples.fasta"
        output_conditioning = f"{output_dir}/{name}-conditioning.pt"
        generate_samples(
            mode='AnalogConditional',
            subset_sequences='-',
            checkpoint_path=args.checkpoint_path,
            output_fasta=output_fasta,
            conditioning_output_path=output_conditioning,
            num_samples=args.num_samples,
            batch_size=args.batch_size,
            length='-',
            charge='-',
            hydrophobicity='-',
            guidance_strength=1.0,
            analog=sequence,
            template='',
            model=model
        )

        # Run prediction on generated samples
        print(f"Running predictions...")
        predict_sequences(
            fasta_file=output_fasta,
            classifier_choice="all",
            output_csv=f"{output_dir}/{name}-predictions.csv",
            predict_proba=True
        )

if __name__ == "__main__":
    main()
