#!/bin/bash

declare -A NEGATIVES=(
    ["GQ20"]="GQLNKFIKKAQRKFHEKFAK"
#    ["Prevotellin-2"]="MLNYLYDRDINRYRAIIKALGLRK"
      ["Bacteroidin-2"]="MGNLVAIVGRPNVGKSTLFNRFH"
#    ["Lachnospirin-1"]="LKQLNRFKYKIIKIKRIIKL"
#    ["Enterococcin-1"]="SIKLIKTVKVVEKIIVEFLKKIKKFSVI"
#    ["Mylodonin-2"]="KRKRGLKLATALSLNNKF"
#    ["Mammuthusin-3"]="KTLKIIRLLF"
    ["Mammutin-1"]="WMTIHALKLSLSFKL"
#    ["Arctodutin-1"]="LLINSIKRLLLGSIL"
#    ["Hesperelin-2"]="KKWVIKSVVFFK"
#    ["Equusin-3"]="KLYQRILWRLISEL"
#    ["Equusin-1"]="FLKLRWSRFARVLL"
    ["DBAASPS_20955"]="IGKLFKRIVERIKRFLRVLLRILR" #hydramp
    ["BoCo1"]="NKIKFINKYVKKVQLKKILVKS" #rampage
    ["DBAASPS_20015"]="LGLFKLLLRLILKGFKL"
    ["DeNo1047"]="ALPSIIKGLLKKL"

)

for NAME in "${!NEGATIVES[@]}"; do
    TEMPLATE="${NEGATIVES[$NAME]}"
    echo $NAME
    echo $TEMPLATE
    OUTPUT_PATH="experiments/results/files/wet-lab/analog-conditioning/timestep-250/analog-subset"
    mkdir -p "${OUTPUT_PATH}"
    CUDA_VISIBLE_DEVICES=3 python project/scripts/inference/generate_samples.py \
        "AdvancedConditional" \
        --subset_sequences data/activity-data/curated-AMPs.fasta \
        --checkpoint_path models/generative_model.ckpt \
        --analog "${TEMPLATE}" \
        --initial_timestep 250 \
        --conditioning_closeness "0.5" \
        --output_fasta "${OUTPUT_PATH}/${NAME}-samples.fasta" \
        --conditioning_output_path "${OUTPUT_PATH}/${NAME}-conditioning.pt" \
        --num_samples 50000 \
        --batch_size 1000 \
        --advanced_conditioning_modes "AnalogConditional,SubsetConditional"

    CUDA_VISIBLE_DEVICES=3 python project/scripts/inference/predict_sequences.py \
        "${OUTPUT_PATH}/${NAME}-samples.fasta" \
        --classifier "all"  \
        --output_csv "${OUTPUT_PATH}/${NAME}-predictions.csv" \
        --get_top_k 500

    python project/scripts/postprocessing/calculate_similarity.py \
        "${TEMPLATE}" \
        "${OUTPUT_PATH}/${NAME}-predictions.csv" \

    python project/scripts/postprocessing/calculate_similarity.py \
        "${TEMPLATE}" \
        "${OUTPUT_PATH}/${NAME}-predictions-best.csv"
done



