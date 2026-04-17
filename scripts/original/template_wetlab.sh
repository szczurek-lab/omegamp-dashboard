#!/bin/bash

declare -A TEMPLATES=(
#    ["cecropine"]="________________KK_____I___I___"
#    ["cecropine_short"]="_____KK_____I___I___"
#    ["sarcotoxin"]="___W_KK________________________________"
#    ["sarcotoxin_short"]="___W_KK_____________"
#    ["p4"]="_________KII__P__K_LL_A__________"
#    ["p4_short"]="___KII__P__K_LL_A___"
    ["bZIP-2"]="___L___V__L___N__L___V__L___V_"
#    ["LPS_Nterminal"]="GWKRKRFG_______"
#    ["LPS_Cterminal"]="_______GWKRKRFG"
)

for NAME in "${!TEMPLATES[@]}"; do
    TEMPLATE="${TEMPLATES[$NAME]}"
    echo $NAME
    echo $TEMPLATE
    OUTPUT_PATH="experiments/results/files/wet-lab/template-conditioning"
    LENGTH=${#TEMPLATE}

#    CUDA_VISIBLE_DEVICES=2 python project/scripts/inference/generate_samples.py \
#        "AdvancedConditional" \
#        --length ${LENGTH} \
#        --charge="2:10" \
#        --hydrophobicity="-0.5:0.8" \
#        --checkpoint_path models/generative_model.ckpt \
#        --template "$TEMPLATE" \
#        --guidance_strength 10 \
#        --output_fasta "${OUTPUT_PATH}/${NAME}-samples.fasta" \
#        --conditioning_output_path "${OUTPUT_PATH}/${NAME}-conditioning.pt" \
#        --num_samples 25000 \
#        --batch_size 1000 \
#        --advanced_conditioning_modes "TemplateConditional,PartialConditional"
#
#
#    python project/scripts/postprocessing/check_template.py \
#      "${OUTPUT_PATH}/${NAME}-samples.fasta" \
#      "${TEMPLATE}" \
#      "${OUTPUT_PATH}/${NAME}-samples-valid.fasta"

    CUDA_VISIBLE_DEVICES=2 python project/scripts/inference/predict_sequences.py \
        "${OUTPUT_PATH}/${NAME}-samples-valid.fasta" \
        --classifier "all"  \
        --output_csv "${OUTPUT_PATH}/${NAME}-predictions.csv" \
        --get_top_k 500

        python project/scripts/postprocessing/calculate_similarity.py \
        "${TEMPLATE}" \
        "${OUTPUT_PATH}/${NAME}-predictions.csv"

done