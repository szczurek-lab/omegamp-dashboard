import os
import argparse
from project.scripts.inference.generate_samples import main as generate_samples
from project.scripts.inference.predict_sequences import main as predict_sequences
from project.config import load_model_for_inference
from project.constants import AMPs_FILE, HQ_AMPs_FILE
from hydra import compose, initialize

TEMPLATES = {
    'cecropine': '________________KK_____I___I___',
    'cecropine_short': '_____KK_____I___I___', # 20 aas
    'sarcotoxin': '___W_KK________________________________',
    'sarcotoxin_short': '___W_KK_____________',
    'p4': '_________KII__P__K_LL_A__________',
    'p4_short': '___KII__P__K_LL_A___',
    'bZIP': '___L___V__L___V__L___V__L___V_',
    'LPS_Nterminal': 'GWKRKRFG_______',
    'LPS_Cterminal': '_______GWKRKRFG',
}

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate sequences using subset conditioning with different subsets')
    parser.add_argument('--num_samples', type=int, default=1000, help='Number of samples to generate for each conditioning')
    parser.add_argument('--batch_size', type=int, default=100, help='Batch size for generation')
    parser.add_argument('--checkpoint_path', type=str, default='models/generative_model.ckpt', help='Path to the model checkpoint')
    parser.add_argument('--output_dir', type=str, default='experiments/results/files/template-conditioning', help='Directory to save output files')

    args = parser.parse_args()

    
    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    with initialize(version_base=None, config_path="../../../config"):
        config = compose(config_name="train")
        model = load_model_for_inference(config, args.checkpoint_path)
        model.to('cuda')

    # Generate samples using regular AMPs file
    for name, template in TEMPLATES.items():
        print(f"Generating {args.num_samples} samples conditioned on template {template}...")
        amp_output = f"{output_dir}/{name}-samples.fasta"
        amp_conditioning = f"{output_dir}/{name}-conditioning.pt"
        generate_samples(
            mode='AdvancedConditional',
            subset_sequences=AMPs_FILE,
            checkpoint_path=args.checkpoint_path,
            output_fasta=amp_output,
            conditioning_output_path=amp_conditioning,
            num_samples=args.num_samples,
            batch_size=args.batch_size,
            length='-',
            charge='-',
            hydrophobicity='-',
            analog='GAAKRAKTAL',
            template=template,
            guidance_strength=1,
            model=model
        )

        # Run prediction on generated samples
        print(f"Running predictions on all samples...")
        predict_sequences(
            fasta_file=amp_output,
            classifier_choice="all",
            output_csv=f"{output_dir}/{name}-predictions.csv",
            predict_proba=True
        )


if __name__ == "__main__":
    main()
