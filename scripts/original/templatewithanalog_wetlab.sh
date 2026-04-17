#!/bin/bash

declare -A TEMPLATES=(
    ["cecropin"]="______________KK_____I___I_____"
    ["sarcotoxin"]="_W_KK__________________________________"
    ["pa4"]="_______KII__P__K_LL_A____________"
    ["bZIP"]="___L___V__L___N__L___V__L___V_"
    ["LG21"]="_____________GWKRKRFG"
)

declare -A PROTOTYPES=(
  ["cecropin"]="SWLSKTAKKLENSAKKRISEGIAIAIQGGPR"
  ["sarcotoxin"]="GWLKKIGKKIERVGQHTRDATIQGLGIAQQAANVAATAR"
  ["pa4"]="GFFALIPKIISSPLFKTLLSAVGSALSSSGGQE"
  ["bZIP"]="MKQLEDKVEELLSKNYHLENEVARLKKLVG"
  ["LG21"]="LLPIVGNLLKSLLGWKRKRFG"
)

CLOSENESS_VALUES=(0 0.25 0.5 0.75 1.0)

for CLOSENESS in "${CLOSENESS_VALUES[@]}"; do
  for NAME in "${!TEMPLATES[@]}"; do
      TEMPLATE="${TEMPLATES[$NAME]}"
      PROTOTYPE="${PROTOTYPES[$NAME]}"
      echo "Processing: $NAME with closeness $CLOSENESS"
      echo $TEMPLATE
      echo $PROTOTYPE
      OUTPUT_PATH="experiments/results/files/wet-lab/templatewithanalog-conditioning/timestep-250/closeness-${CLOSENESS}"
      mkdir -p "${OUTPUT_PATH}"

      LENGTH=${#TEMPLATE}

      CUDA_VISIBLE_DEVICES=2 python project/scripts/inference/generate_samples.py \
          "AdvancedConditional" \
          --length ${LENGTH} \
          --checkpoint_path models/generative_model.ckpt \
          --analog "${PROTOTYPE}" \
          --template "${TEMPLATE}" \
          --initial_timestep 250 \
          --conditioning_closeness "${CLOSENESS}"  \
          --guidance_strength 1 \
          --output_fasta "${OUTPUT_PATH}/${NAME}-samples.fasta" \
          --conditioning_output_path "${OUTPUT_PATH}/${NAME}-conditioning.pt" \
          --num_samples 50000 \
          --batch_size 1000 \
          --advanced_conditioning_modes "TemplateConditional,AnalogConditional,PartialConditional"

      python project/scripts/postprocessing/check_template.py \
        "${OUTPUT_PATH}/${NAME}-samples.fasta" \
        "${TEMPLATE}" \
        "${OUTPUT_PATH}/${NAME}-samples-valid.fasta"

      CUDA_VISIBLE_DEVICES=2 python project/scripts/inference/predict_sequences.py \
          "${OUTPUT_PATH}/${NAME}-samples-valid.fasta" \
          --classifier "all"  \
          --output_csv "${OUTPUT_PATH}/${NAME}-predictions.csv" \
          --get_top_k 500

      python project/scripts/postprocessing/calculate_similarity.py \
          "${PROTOTYPE}" \
          "${OUTPUT_PATH}/${NAME}-predictions.csv"

      done
  done