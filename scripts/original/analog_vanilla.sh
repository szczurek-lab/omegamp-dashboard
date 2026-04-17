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

CLOSENESS_VALUES=(0 0.25 0.5 0.75 1)

for CLOSENESS in "${CLOSENESS_VALUES[@]}"; do
    for NAME in "${!NEGATIVES[@]}"; do
        ANALOG="${NEGATIVES[$NAME]}"
        echo "Processing: $NAME with closeness $CLOSENESS"
        echo "Analog: $ANALOG"

        OUTPUT_PATH="experiments/results/files/wet-lab/analog-conditioning/timestep-100/closeness-${CLOSENESS}"

        mkdir -p "${OUTPUT_PATH}"

        CUDA_VISIBLE_DEVICES=1 python project/scripts/inference/generate_samples.py \
            "AnalogConditional" \
            --checkpoint_path models/generative_model.ckpt \
            --analog "${ANALOG}" \
            --conditioning_closeness "${CLOSENESS}" \
            --initial_timestep 100 \
            --output_fasta "${OUTPUT_PATH}/${NAME}-samples.fasta" \
            --conditioning_output_path "${OUTPUT_PATH}/${NAME}-conditioning.pt" \
            --num_samples 50000 \
            --batch_size 1000 \
            --advanced_conditioning_modes "AnalogConditional"

        CUDA_VISIBLE_DEVICES=1 python project/scripts/inference/predict_sequences.py \
            "${OUTPUT_PATH}/${NAME}-samples.fasta" \
            --classifier "all"  \
            --output_csv "${OUTPUT_PATH}/${NAME}-predictions.csv" \
            --get_top_k 500

    done
done