import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import pandas as pd
from io import StringIO

def extract_blocks_from_alignment(
    alignment,
    min_column_completeness=0.9,
    min_block_length=2
):
    """
    Extract ungapped alignment blocks from a BioPython MultipleSeqAlignment object.
    """
    aln_array = pd.DataFrame(
        [list(str(rec.seq)) for rec in alignment],
        index=[rec.id for rec in alignment]
    )

    # Identify valid columns
    non_gap_fraction = (aln_array != "-").sum(axis=0) / aln_array.shape[0]
    valid_columns = non_gap_fraction >= min_column_completeness

    # Group valid columns into blocks
    blocks = []
    current = []
    for i, is_valid in enumerate(valid_columns):
        if is_valid:
            current.append(i)
        elif current:
            if len(current) >= min_block_length:
                blocks.append(current)
            current = []
    if current and len(current) >= min_block_length:
        blocks.append(current)

    # Extract valid blocks
    msa_blocks = []
    for cols in blocks:
        block_df = aln_array.iloc[:, cols]
        no_gap_rows = block_df[(block_df != "-").all(axis=1)]

        if len(no_gap_rows) >= 2:
            records = [
                SeqRecord(Seq("".join(row)), id=row_id, description="")
                for row_id, row in no_gap_rows.iterrows()
            ]
            block_alignment = MultipleSeqAlignment(records)
            msa_blocks.append(block_alignment)

    return msa_blocks


def extract_blocks_from_fasta_string(
    fasta_str,
    min_column_completeness=0.9,
    min_block_length=2
):
    """
    Extract blocks from a FASTA alignment string (in-memory).
    """
    alignment = AlignIO.read(StringIO(fasta_str), "fasta")
    blocks = extract_blocks_from_alignment(alignment, min_column_completeness, min_block_length)

    # Convert each block back to FASTA string
    fasta_blocks = []
    for block in blocks:
        output = StringIO()
        AlignIO.write(block, output, "fasta")
        fasta_blocks.append(output.getvalue())
    return fasta_blocks


def extract_blocks_to_files(
    input_path,
    output_dir,
    min_column_completeness=0.9,
    min_block_length=2
):
    """
    Original file-based implementation. Saves blocks as separate FASTA files.
    """
    os.makedirs(output_dir, exist_ok=True)
    alignment = AlignIO.read(input_path, "fasta")
    blocks = extract_blocks_from_alignment(alignment, min_column_completeness, min_block_length)

    for i, block in enumerate(blocks):
        output_file = os.path.join(output_dir, f"block_{i+1}.fasta")
        AlignIO.write(block, output_file, "fasta")
        print(f"[INFO] Saved block {i+1}: {output_file}")
