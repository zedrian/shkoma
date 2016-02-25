def find(protein_sequence, peptide_sequence):
    peptide_without_ends_sequence = peptide_sequence[1:-2]
    start = protein_sequence.find(peptide_without_ends_sequence)

    if start == -1:
        return -1

    return start - 1
