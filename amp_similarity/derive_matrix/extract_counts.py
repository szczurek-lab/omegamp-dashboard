from pysum.src.algorithm import blosum

def load_sequences(filename):
    with open(filename, 'r') as f:
        sequences = [line.strip() for line in f if line.strip()]
    return sequences

def run_blosum_from_seqs(plain_str: str, matrix_version: int = 80) -> dict:
    
    seq_array = [line.strip() for line in plain_str.strip().splitlines() if line.strip()]

    counts = blosum(seq_array, matrix_version)

    if counts is None:
        print("Only one sequence found, no counts to return.")
        return {}

    return counts