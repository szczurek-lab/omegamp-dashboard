import subprocess

def run_mafft(fasta_str, matrix_path="custom_aa_matrix.txt"):
    """
    Run MAFFT with a custom substitution matrix on an in-memory FASTA string.

    Returns:
        str: Aligned FASTA as a string.
    """
    try:
        result = subprocess.run(
            [
                "mafft",
                "--amino",
                "--textmatrix", matrix_path, # Comment if using a blossum62
                "--localpair",
                "-"
            ],
            input=fasta_str.encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=True
        )
        return result.stdout.decode()
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] MAFFT failed: {e}")
        return ""
