#!/bin/bash

# Biologically-motivated motif sampling
# Template length: 25 residues
# Positions: N-term (0-4), Central (10-14), C-term (20-24)

OUTPUT_PATH="experiments/results/files/template_sweep"
mkdir -p ${OUTPUT_PATH}

CHECKPOINT="models/generative_model.ckpt"
NUM_SAMPLES=2000
BATCH_SIZE=500
GUIDANCE=10

# Function to run generation
run_generation() {
    local NAME=$1
    local TEMPLATE=$2
    local LENGTH=${#TEMPLATE}

    echo "========================================="
    echo "Generating: ${NAME}"
    echo "Template: ${TEMPLATE}"
    echo "Length: ${LENGTH}"
    echo "========================================="

    python project/scripts/inference/generate_samples.py \
        "AdvancedConditional" \
        --length ${LENGTH} \
        --charge="2:10" \
        --hydrophobicity="-0.5:0.8" \
        --checkpoint_path ${CHECKPOINT} \
        --template "${TEMPLATE}" \
        --guidance_strength ${GUIDANCE} \
        --output_fasta "${OUTPUT_PATH}/${NAME}.fasta" \
        --conditioning_output_path "${OUTPUT_PATH}/${NAME}_conditioning.pt" \
        --num_samples ${NUM_SAMPLES} \
        --batch_size ${BATCH_SIZE} \
        --advanced_conditioning_modes "TemplateConditional,PartialConditional"
}

################################################################################
# CONTIGUOUS MOTIFS - CHARGED (Short: 3-4, Medium: 5-6)
################################################################################

## Short (3-4 residues) - Charged
# N-terminal
run_generation "cont_charged_short_KKK_Nterm" "KKK______________________"
run_generation "cont_charged_short_KKKR_Nterm" "KKKR_____________________"
run_generation "cont_charged_short_RRRK_Nterm" "RRRK_____________________"

# Central
run_generation "cont_charged_short_KKK_Central" "___________KKK___________"
run_generation "cont_charged_short_KKKR_Central" "___________KKKR__________"
run_generation "cont_charged_short_RRRK_Central" "___________RRRK__________"

# C-terminal
run_generation "cont_charged_short_KKK_Cterm" "______________________KKK"
run_generation "cont_charged_short_KKKR_Cterm" "_____________________KKKR"
run_generation "cont_charged_short_RRRK_Cterm" "_____________________RRRK"

## Medium (5-6 residues) - Charged
# N-terminal
run_generation "cont_charged_med_KKKKK_Nterm" "KKKKK____________________"
run_generation "cont_charged_med_KRKRKR_Nterm" "KRKRKR___________________"

# Central
run_generation "cont_charged_med_KKKKK_Central" "__________KKKKK__________"
run_generation "cont_charged_med_KRKRKR_Central" "__________KRKRKR_________"

# C-terminal
run_generation "cont_charged_med_KKKKK_Cterm" "____________________KKKKK"
run_generation "cont_charged_med_KRKRKR_Cterm" "___________________KRKRKR"

################################################################################
# CONTIGUOUS MOTIFS - HYDROPHOBIC (Short: 3-4, Medium: 5-6)
################################################################################

## Short (3-4 residues) - Hydrophobic
# N-terminal
run_generation "cont_hydrophobic_short_FFF_Nterm" "FFF______________________"
run_generation "cont_hydrophobic_short_LLLL_Nterm" "LLLL_____________________"
run_generation "cont_hydrophobic_short_IIII_Nterm" "IIII_____________________"

# Central
run_generation "cont_hydrophobic_short_FFF_Central" "___________FFF___________"
run_generation "cont_hydrophobic_short_LLLL_Central" "___________LLLL__________"
run_generation "cont_hydrophobic_short_IIII_Central" "___________IIII__________"

# C-terminal
run_generation "cont_hydrophobic_short_FFF_Cterm" "______________________FFF"
run_generation "cont_hydrophobic_short_LLLL_Cterm" "_____________________LLLL"
run_generation "cont_hydrophobic_short_IIII_Cterm" "_____________________IIII"

## Medium (5-6 residues) - Hydrophobic
# N-terminal
run_generation "cont_hydrophobic_med_FFFFF_Nterm" "FFFFF____________________"
run_generation "cont_hydrophobic_med_LLLLLL_Nterm" "LLLLLL___________________"

# Central
run_generation "cont_hydrophobic_med_FFFFF_Central" "__________FFFFF__________"
run_generation "cont_hydrophobic_med_LLLLLL_Central" "__________LLLLLL_________"

# C-terminal
run_generation "cont_hydrophobic_med_FFFFF_Cterm" "____________________FFFFF"
run_generation "cont_hydrophobic_med_LLLLLL_Cterm" "___________________LLLLLL"

################################################################################
# CONTIGUOUS MOTIFS - AMPHIPATHIC (Short: 4, Medium: 6)
################################################################################

## Short (4 residues) - Amphipathic
# N-terminal
run_generation "cont_amphipathic_short_KFKF_Nterm" "KFKF_____________________"
run_generation "cont_amphipathic_short_KWKW_Nterm" "KWKW_____________________"
run_generation "cont_amphipathic_short_KLKL_Nterm" "KLKL_____________________"

# Central
run_generation "cont_amphipathic_short_KFKF_Central" "___________KFKF__________"
run_generation "cont_amphipathic_short_KWKW_Central" "___________KWKW__________"
run_generation "cont_amphipathic_short_KLKL_Central" "___________KLKL__________"

# C-terminal
run_generation "cont_amphipathic_short_KFKF_Cterm" "_____________________KFKF"
run_generation "cont_amphipathic_short_KWKW_Cterm" "_____________________KWKW"
run_generation "cont_amphipathic_short_KLKL_Cterm" "_____________________KLKL"

## Medium (6 residues) - Amphipathic
# N-terminal
run_generation "cont_amphipathic_med_KFKFKF_Nterm" "KFKFKF___________________"
run_generation "cont_amphipathic_med_KLKLKL_Nterm" "KLKLKL___________________"

# Central
run_generation "cont_amphipathic_med_KFKFKF_Central" "__________KFKFKF_________"
run_generation "cont_amphipathic_med_KLKLKL_Central" "__________KLKLKL_________"

# C-terminal
run_generation "cont_amphipathic_med_KFKFKF_Cterm" "___________________KFKFKF"
run_generation "cont_amphipathic_med_KLKLKL_Cterm" "___________________KLKLKL"

################################################################################
# CONTIGUOUS MOTIFS - STRUCTURAL
################################################################################

## Glycine-rich (flexibility) - Short
# N-terminal
run_generation "cont_structural_short_GG_Nterm" "GG_______________________"
run_generation "cont_structural_short_GGG_Nterm" "GGG______________________"

# Central
run_generation "cont_structural_short_GG_Central" "____________GG___________"
run_generation "cont_structural_short_GGG_Central" "___________GGG___________"

# C-terminal
run_generation "cont_structural_short_GG_Cterm" "_______________________GG"
run_generation "cont_structural_short_GGG_Cterm" "______________________GGG"

## Proline-rich (turns) - Short only
# N-terminal
run_generation "cont_structural_short_PP_Nterm" "PP_______________________"
run_generation "cont_structural_short_PPP_Nterm" "PPP______________________"

# Central
run_generation "cont_structural_short_PP_Central" "____________PP___________"
run_generation "cont_structural_short_PPP_Central" "___________PPP___________"

# C-terminal
run_generation "cont_structural_short_PP_Cterm" "_______________________PP"
run_generation "cont_structural_short_PPP_Cterm" "______________________PPP"

## Neutral/Polar - Short (using S, T, A)
# N-terminal
run_generation "cont_neutral_short_SSSS_Nterm" "SSSS_____________________"
run_generation "cont_neutral_short_AAAA_Nterm" "AAAA_____________________"

# Central
run_generation "cont_neutral_short_SSSS_Central" "___________SSSS__________"
run_generation "cont_neutral_short_AAAA_Central" "___________AAAA__________"

# C-terminal
run_generation "cont_neutral_short_SSSS_Cterm" "_____________________SSSS"
run_generation "cont_neutral_short_AAAA_Cterm" "_____________________AAAA"

################################################################################
# SPACED MOTIFS - SHORT SPAN
################################################################################

## Charged - Spaced (span 5-7)
# N-terminal
run_generation "spaced_charged_K_K_Nterm" "K___K____________________"
run_generation "spaced_charged_KK_KK_Nterm" "KK__KK___________________"
run_generation "spaced_charged_K_R_K_Nterm" "K__R__K__________________"

# Central
run_generation "spaced_charged_K_K_Central" "__________K___K__________"
run_generation "spaced_charged_KK_KK_Central" "_________KK__KK__________"
run_generation "spaced_charged_K_R_K_Central" "_________K__R__K_________"

# C-terminal
run_generation "spaced_charged_K_K_Cterm" "____________________K___K"
run_generation "spaced_charged_KK_KK_Cterm" "___________________KK__KK"
run_generation "spaced_charged_K_R_K_Cterm" "__________________K__R__K"

## Hydrophobic - Spaced (span 5-7)
# N-terminal
run_generation "spaced_hydrophobic_F_F_Nterm" "F___F____________________"
run_generation "spaced_hydrophobic_L_L_L_Nterm" "L__L__L__________________"
run_generation "spaced_hydrophobic_W_W_Nterm" "W____W___________________"

# Central
run_generation "spaced_hydrophobic_F_F_Central" "__________F___F__________"
run_generation "spaced_hydrophobic_L_L_L_Central" "_________L__L__L_________"
run_generation "spaced_hydrophobic_W_W_Central" "_________W____W__________"

# C-terminal
run_generation "spaced_hydrophobic_F_F_Cterm" "____________________F___F"
run_generation "spaced_hydrophobic_L_L_L_Cterm" "__________________L__L__L"
run_generation "spaced_hydrophobic_W_W_Cterm" "___________________W____W"

## Disulfide bridges - Spaced
# N-terminal
run_generation "spaced_structural_C_C_short_Nterm" "C____C___________________"
run_generation "spaced_structural_C_C_long_Nterm" "C________C_______________"

# Central
run_generation "spaced_structural_C_C_short_Central" "_________C____C__________"
run_generation "spaced_structural_C_C_long_Central" "________C________C_______"

# C-terminal
run_generation "spaced_structural_C_C_short_Cterm" "___________________C____C"
run_generation "spaced_structural_C_C_long_Cterm" "_______________C________C"

################################################################################
# SPACED MOTIFS - MEDIUM/LONG SPAN
################################################################################

## Charged - Long span (span 9-13)
# N-terminal
run_generation "spaced_charged_K_K_longspan_Nterm" "K________K_______________"
run_generation "spaced_charged_K_K_K_longspan_Nterm" "K____K____K______________"

# Central
run_generation "spaced_charged_K_K_longspan_Central" "________K________K_______"
run_generation "spaced_charged_K_K_K_longspan_Central" "_______K____K____K_______"

# C-terminal
run_generation "spaced_charged_K_K_longspan_Cterm" "_______________K________K"
run_generation "spaced_charged_K_K_K_longspan_Cterm" "______________K____K____K"

## Hydrophobic - Long span (span 9-13)
# N-terminal
run_generation "spaced_hydrophobic_F_F_longspan_Nterm" "F________F_______________"
run_generation "spaced_hydrophobic_L_L_L_longspan_Nterm" "L____L____L______________"

# Central
run_generation "spaced_hydrophobic_F_F_longspan_Central" "________F________F_______"
run_generation "spaced_hydrophobic_L_L_L_longspan_Central" "_______L____L____L_______"

# C-terminal
run_generation "spaced_hydrophobic_F_F_longspan_Cterm" "_______________F________F"
run_generation "spaced_hydrophobic_L_L_L_longspan_Cterm" "______________L____L____L"

################################################################################
# BIOLOGICALLY-KNOWN MOTIFS
################################################################################

## LPS/Heparin-binding motifs (basic-rich clusters)
# Classic LPS-binding: GWKRKRFG
run_generation "bio_LPS_binding_Nterm" "GWKRKRFG_________________"
run_generation "bio_LPS_binding_Central" "_________GWKRKRFG________"
run_generation "bio_LPS_binding_Cterm" "_________________GWKRKRFG"

# Simplified LPS-binding variants
run_generation "bio_LPS_simple1_Nterm" "KWKRKR___________________"
run_generation "bio_LPS_simple1_Central" "__________KWKRKR_________"
run_generation "bio_LPS_simple2_Nterm" "GKRKR____________________"
run_generation "bio_LPS_simple2_Central" "___________GKRKR_________"

## Cell-penetrating peptide-like motifs (basic clusters)
# TAT-like: RKKRRQRRR (simplified to avoid Q)
run_generation "bio_CPP_like1_Nterm" "RKKRRR___________________"
run_generation "bio_CPP_like1_Central" "__________RKKRRR_________"
run_generation "bio_CPP_like2_Nterm" "KRRKRR___________________"

## α-helical antimicrobial patterns
# KLAL repeat (magainin-like)
run_generation "bio_alpha_helix_KLAL_Nterm" "KLALKLAL_________________"
run_generation "bio_alpha_helix_KLAL_Central" "_________KLALKLAL________"
run_generation "bio_alpha_helix_KLAL_Cterm" "_________________KLALKLAL"

# KWKS repeat
run_generation "bio_alpha_helix_KWKS_Nterm" "KWKSKWKS_________________"
run_generation "bio_alpha_helix_KWKS_Central" "_________KWKSKWKS________"

## Membrane-disrupting patterns (alternating charged/hydrophobic)
run_generation "bio_membrane_KLFKLF_Nterm" "KLFKLF___________________"
run_generation "bio_membrane_KLFKLF_Central" "__________KLFKLF_________"
run_generation "bio_membrane_KIFKIF_Nterm" "KIFKIF___________________"
run_generation "bio_membrane_KIFKIF_Central" "__________KIFKIF_________"

## Metal-binding motifs (His-rich)
# HX3H pattern
run_generation "bio_metal_HxxxH_Nterm" "H___H____________________"
run_generation "bio_metal_HxxxH_Central" "__________H___H__________"
run_generation "bio_metal_HxxxH_Cterm" "____________________H___H"

# HX4H pattern
run_generation "bio_metal_HxxxxH_Nterm" "H____H___________________"
run_generation "bio_metal_HxxxxH_Central" "_________H____H__________"

## β-sheet antimicrobial patterns (defensin-like)
# Alternating charged/hydrophobic with spacing
run_generation "bio_beta_sheet_pattern1_Nterm" "K_F_K_F__________________"
run_generation "bio_beta_sheet_pattern1_Central" "_________K_F_K_F_________"
run_generation "bio_beta_sheet_pattern2_Nterm" "R_L_R_L__________________"

## Tryptophan-rich antimicrobial motifs
run_generation "bio_Trp_rich_WKW_Nterm" "WKW______________________"
run_generation "bio_Trp_rich_WKW_Central" "___________WKW___________"
run_generation "bio_Trp_rich_WKWK_Nterm" "WKWK_____________________"
run_generation "bio_Trp_rich_WKWKW_Nterm" "WKWKW____________________"

## Cecropin-like patterns
run_generation "bio_cecropin_like1_Nterm" "KWKLFKKI_________________"
run_generation "bio_cecropin_like1_Central" "_________KWKLFKKI________"

## Cathelicidin-like patterns (LLGDFFR-like)
run_generation "bio_cathelicidin_like_Nterm" "LLGDFFR__________________"
run_generation "bio_cathelicidin_like_Central" "_________LLGDFFR_________"

################################################################################
# VARIABLE POSITION (test model flexibility across lengths)
################################################################################

## Test same motif at different sequence lengths (15, 20, 30)
# Short sequences (15 aa)
run_generation "var_length15_KKK_pos3" "__KKK__________"
run_generation "var_length15_KWKW_pos5" "____KWKW_______"

# Medium sequences (20 aa)
run_generation "var_length20_KKK_pos3" "__KKK_______________"
run_generation "var_length20_KWKW_pos8" "_______KWKW_________"

# Long sequences (30 aa)
run_generation "var_length30_KKK_pos3" "__KKK_________________________"
run_generation "var_length30_KWKW_pos13" "____________KWKW______________"

echo "========================================="
echo "Template sweep complete!"
echo "Total motifs tested: ~130"
echo "Results saved to: ${OUTPUT_PATH}"
echo "========================================="