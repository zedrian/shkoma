def find(protein_sequence, peptide_sequence):
    peptide_without_ends_sequence = peptide_sequence[1:-2]
    start = protein_sequence.find(peptide_without_ends_sequence)

    if start == -1:
        return -1

    return start - 1


# returns list of sequences remaining after cutting peptide sequences from protein one
def cut_received_peptide_sequences(protein_sequence, peptide_sequences):
    modified_protein_sequence = protein_sequence

    for peptide_sequence in peptide_sequences:
        start = find(modified_protein_sequence, peptide_sequence)
        if start != -1:
            line = ''
            for i in range(0, len(modified_protein_sequence)):
                if i not in range(start, start + len(peptide_sequence)):
                    line += modified_protein_sequence[i]
                else:
                    line += '-'
            modified_protein_sequence = line

    missed_peptide_sequences = modified_protein_sequence.split('-')
    missed_peptide_sequences = [fragment for fragment in missed_peptide_sequences if fragment != '']
    return missed_peptide_sequences
