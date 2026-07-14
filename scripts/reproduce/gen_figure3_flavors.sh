#!/usr/bin/env bash
# Figure 3 flavors: nine OmegAMP generation modes used in the figure 3
# comparison. Score separately with score_figure3_flavors.sh.
# Run under the OmegAMP env.
# Usage: bash scripts/reproduce/gen_figure3_flavors.sh [--mini] [--only positive|negative]
set -euo pipefail

source config.sh

DENOVO_SAMPLES=50000
ANALOGS_PER_PROTO=10
LENGTH="5:30"
CHARGE="2:10"
HYDRO="-0.5:0.8"
MOTIF="G----G"
MOTIF_SEED=42

MINI_PROTOS=0
ONLY=both
while [ $# -gt 0 ]; do
    case "$1" in
        --mini)
            DENOVO_SAMPLES=100
            ANALOGS_PER_PROTO=3
            MINI_PROTOS=5
            ;;
        --only) ONLY="$2"; shift ;;
    esac
    shift
done

OUT="${DATA_DIR}/figure3/flavors"
WORK="${OUT}/_inputs"
mkdir -p "${OUT}" "${WORK}"

# In mini mode, subsample to first MINI_PROTOS prototypes for the analog runs.
prototype_subset() {
    local src="$1" out="$2" n="$3"
    awk -v n="${n}" '/^>/{c++; if (c > n) exit} {print}' "${src}" > "${out}"
}

# Insert MOTIF at a seeded random position in each prototype, producing one
# length-matched dashed template per prototype, in prototype order.
make_motif_templates() {
    local proto_fasta="$1" out_fasta="$2"
    python - "${proto_fasta}" "${out_fasta}" "${MOTIF}" "${MOTIF_SEED}" <<'PY'
import sys, random
proto, out, motif, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
rng = random.Random(seed)
entries, hdr, seq = [], None, []
for line in open(proto):
    line = line.rstrip()
    if line.startswith(">"):
        if hdr is not None:
            entries.append((hdr, "".join(seq)))
        hdr, seq = line, []
    else:
        seq.append(line)
if hdr is not None:
    entries.append((hdr, "".join(seq)))
with open(out, "w") as f:
    for hdr, s in entries:
        L = len(s)
        if L < len(motif):
            t = motif[:L]
        else:
            start = rng.randint(0, L - len(motif))
            t = "-" * start + motif + "-" * (L - start - len(motif))
        f.write(f"{hdr}\n{t}\n")
PY
}

cd "${OMEGAMP_DIR}"

# De novo flavors. No prototype set, so --only doesn't apply.

fasta="${OUT}/omegamp_du/OmegAMP-U.fasta"
mkdir -p "$(dirname "${fasta}")"
if [ ! -s "${fasta}" ]; then
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo unconditional \
        --checkpoint_path "${CHECKPOINT}" \
        --num_samples "${DENOVO_SAMPLES}" --batch_size 256 \
        --output_fasta "${fasta}" \
        --conditioning_output_path "${OUT}/omegamp_du/conditioning.pt"
fi

fasta="${OUT}/omegamp_dt/OmegAMP-T.fasta"
mkdir -p "$(dirname "${fasta}")"
if [ ! -s "${fasta}" ]; then
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
        --checkpoint_path "${CHECKPOINT}" \
        --length="${LENGTH}" --charge="${CHARGE}" --hydrophobicity="${HYDRO}" \
        --num_samples "${DENOVO_SAMPLES}" --batch_size 256 \
        --output_fasta "${fasta}" \
        --conditioning_output_path "${OUT}/omegamp_dt/conditioning.pt"
fi

fasta="${OUT}/omegamp_dp/denovo_subset.fasta"
mkdir -p "$(dirname "${fasta}")"
if [ ! -s "${fasta}" ]; then
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo prototype-derived \
        --checkpoint_path "${CHECKPOINT}" \
        --prototype_sequences "${POSITIVE_PROTOTYPES}" \
        --sigma "${OPTIMAL_POSITIVE_SIGMA}" \
        --num_samples "${DENOVO_SAMPLES}" --batch_size 256 \
        --output_fasta "${fasta}" \
        --conditioning_output_path "${OUT}/omegamp_dp/conditioning.pt"
fi

# Prototype baselines: copy the source FASTAs.
mkdir -p "${OUT}/prototypes"
[ "${ONLY}" != "negative" ] && cp "${POSITIVE_PROTOTYPES}" "${OUT}/prototypes/positive_prototypes.fasta"
[ "${ONLY}" != "positive" ] && cp "${NEGATIVE_PROTOTYPES}" "${OUT}/prototypes/negative_prototypes.fasta"

# Analog flavors. Five modes per prototype set, sharing two repeated-input
# FASTAs (one prototypes-repeated, one motif-templates-repeated).
run_analog_set() {
    local pt="$1" proto_fasta="$2" tau="$3" sigma="$4"

    if [ "${MINI_PROTOS}" -gt 0 ]; then
        local mini_proto="${WORK}/${pt}_mini_protos.fasta"
        prototype_subset "${proto_fasta}" "${mini_proto}" "${MINI_PROTOS}"
        proto_fasta="${mini_proto}"
    fi

    local proto_rep="${WORK}/${pt}_proto_repeated.fasta"
    local motif_rep="${WORK}/${pt}_motif_repeated.fasta"
    [ -s "${proto_rep}" ] || repeat_fasta "${proto_fasta}" "${ANALOGS_PER_PROTO}" "${proto_rep}"
    if [ ! -s "${motif_rep}" ]; then
        make_motif_templates "${proto_fasta}" "${WORK}/${pt}_motif_templates.fasta"
        repeat_fasta "${WORK}/${pt}_motif_templates.fasta" "${ANALOGS_PER_PROTO}" "${motif_rep}"
    fi
    local n=$(grep -c '^>' "${proto_rep}")

    local fasta="${OUT}/omegamp_au/${pt}_analog_only.fasta"
    mkdir -p "$(dirname "${fasta}")"
    if [ ! -s "${fasta}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog unconditional \
            --checkpoint_path "${CHECKPOINT}" \
            --analog_sequences "${proto_rep}" \
            --tau "${tau}" \
            --num_samples "${n}" --batch_size 256 \
            --output_fasta "${fasta}" \
            --conditioning_output_path "${OUT}/omegamp_au/${pt}_conditioning.pt"
    fi

    fasta="${OUT}/omegamp_at/${pt}_analog_property.fasta"
    mkdir -p "$(dirname "${fasta}")"
    if [ ! -s "${fasta}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog targeted \
            --checkpoint_path "${CHECKPOINT}" \
            --analog_sequences "${proto_rep}" \
            --tau "${tau}" \
            --length="${LENGTH}" --charge="${CHARGE}" --hydrophobicity="${HYDRO}" \
            --num_samples "${n}" --batch_size 256 \
            --output_fasta "${fasta}" \
            --conditioning_output_path "${OUT}/omegamp_at/${pt}_conditioning.pt"
    fi

    fasta="${OUT}/omegamp_am/${pt}_analog_template.fasta"
    mkdir -p "$(dirname "${fasta}")"
    if [ ! -s "${fasta}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog-motif unconditional \
            --checkpoint_path "${CHECKPOINT}" \
            --analog_sequences "${proto_rep}" \
            --motif_sequences "${motif_rep}" \
            --tau "${tau}" --guidance_strength 1.0 \
            --num_samples "${n}" --batch_size 256 \
            --output_fasta "${fasta}" \
            --conditioning_output_path "${OUT}/omegamp_am/${pt}_conditioning.pt"
    fi

    fasta="${OUT}/omegamp_amt/${pt}_analog_template_property.fasta"
    mkdir -p "$(dirname "${fasta}")"
    if [ ! -s "${fasta}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog-motif targeted \
            --checkpoint_path "${CHECKPOINT}" \
            --analog_sequences "${proto_rep}" \
            --motif_sequences "${motif_rep}" \
            --tau "${tau}" --guidance_strength 1.0 \
            --length="${LENGTH}" --charge="${CHARGE}" --hydrophobicity="${HYDRO}" \
            --num_samples "${n}" --batch_size 256 \
            --output_fasta "${fasta}" \
            --conditioning_output_path "${OUT}/omegamp_amt/${pt}_conditioning.pt"
    fi

    fasta="${OUT}/omegamp_ap/${pt}_analog_subset.fasta"
    mkdir -p "$(dirname "${fasta}")"
    if [ ! -s "${fasta}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog prototype-derived \
            --checkpoint_path "${CHECKPOINT}" \
            --analog_sequences "${proto_rep}" \
            --prototype_sequences "${proto_rep}" \
            --tau "${tau}" --sigma "${sigma}" \
            --num_samples "${n}" --batch_size 256 \
            --output_fasta "${fasta}" \
            --conditioning_output_path "${OUT}/omegamp_ap/${pt}_conditioning.pt"
    fi
}

[ "${ONLY}" != "negative" ] && run_analog_set positive "${POSITIVE_PROTOTYPES}" "${OPTIMAL_POSITIVE_TAU}" "${OPTIMAL_POSITIVE_SIGMA}"
[ "${ONLY}" != "positive" ] && run_analog_set negative "${NEGATIVE_PROTOTYPES}" "${OPTIMAL_NEGATIVE_TAU}" "${OPTIMAL_NEGATIVE_SIGMA}"

echo "Done. $(find "${OUT}" -name '*.fasta' -not -path '*/_inputs/*' | wc -l) FASTA files in ${OUT}"