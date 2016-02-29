from Bio.SeqUtils.ProtParam import ProteinAnalysis


class ProteinParameters:
    def __init__(self, sequence):
        self.sequence_length = len(sequence)
        analysis = ProteinAnalysis(sequence)

        self.amino_acids_percents = analysis.get_amino_acids_percent()
        self.amino_acids_composition = calculate_amino_acids_composition(sequence)
        self.aromaticity = analysis.aromaticity()
        self.instability = analysis.instability_index()
        self.flexibility = analysis.flexibility()
        self.weight_list = None  # analysis.weight_list(11, ?)
        self.protein_scale = None  # analysis.protein_scale(?, 11, ?)
        self.isoelectric_point = analysis.isoelectric_point()
        self.secondary_structure_fraction = analysis.secondary_structure_fraction()
        self.molecular_weight = analysis.molecular_weight()
        self.kyte_plot = analysis.gravy()
        self.M = 0
        self.Z = 0

    def __str__(self):
        return '  Sequence length: ' + str(self.sequence_length) + '\n' + \
               '  Amino acids percents:\n' \
               '    ACID\t\tPERCENT\n' + \
               '    A\t\t\t' + str(self.amino_acids_percents['A'] * 100.0) + '\n' + \
               '    G\t\t\t' + str(self.amino_acids_percents['G'] * 100.0) + '\n' + \
               '    V\t\t\t' + str(self.amino_acids_percents['V'] * 100.0) + '\n' + \
               '    M\t\t\t' + str(self.amino_acids_percents['M'] * 100.0) + '\n' + \
               '    D\t\t\t' + str(self.amino_acids_percents['D'] * 100.0) + '\n' + \
               '    Y\t\t\t' + str(self.amino_acids_percents['Y'] * 100.0) + '\n' + \
               '    N\t\t\t' + str(self.amino_acids_percents['N'] * 100.0) + '\n' + \
               '    S\t\t\t' + str(self.amino_acids_percents['S'] * 100.0) + '\n' + \
               '    W\t\t\t' + str(self.amino_acids_percents['W'] * 100.0) + '\n' + \
               '    L\t\t\t' + str(self.amino_acids_percents['L'] * 100.0) + '\n' + \
               '    F\t\t\t' + str(self.amino_acids_percents['F'] * 100.0) + '\n' + \
               '    I\t\t\t' + str(self.amino_acids_percents['I'] * 100.0) + '\n' + \
               '    K\t\t\t' + str(self.amino_acids_percents['K'] * 100.0) + '\n' + \
               '    P\t\t\t' + str(self.amino_acids_percents['P'] * 100.0) + '\n' + \
               '    Q\t\t\t' + str(self.amino_acids_percents['Q'] * 100.0) + '\n' + \
               '    C\t\t\t' + str(self.amino_acids_percents['C'] * 100.0) + '\n' + \
               '    E\t\t\t' + str(self.amino_acids_percents['E'] * 100.0) + '\n' + \
               '    R\t\t\t' + str(self.amino_acids_percents['R'] * 100.0) + '\n' + \
               '    T\t\t\t' + str(self.amino_acids_percents['T'] * 100.0) + '\n' + \
               '    H\t\t\t' + str(self.amino_acids_percents['H'] * 100.0) + '\n' + \
               '  Amino acids composition:\n' + \
               '    GROUP\t\tNUMBER\t\tPERCENT\n' + \
               '    Small\t\t' + str(self.amino_acids_composition['Small']) + '\t\t\t' + str(
            self.amino_acids_composition['Small'] / self.sequence_length * 100.0) + '\n' + \
               '    Aliphatic\t' + str(self.amino_acids_composition['Aliphatic']) + '\t\t\t' + str(
            self.amino_acids_composition['Aliphatic'] / self.sequence_length * 100.0) + '\n' + \
               '    Aromatic\t' + str(self.amino_acids_composition['Aromatic']) + '\t\t\t' + str(
            self.amino_acids_composition['Aromatic'] / self.sequence_length * 100.0) + '\n' + \
               '    Non-polar\t' + str(self.amino_acids_composition['Non-polar']) + '\t\t\t' + str(
            self.amino_acids_composition['Non-polar'] / self.sequence_length * 100.0) + '\n' + \
               '    Polar\t\t' + str(self.amino_acids_composition['Polar']) + '\t\t\t' + str(
            self.amino_acids_composition['Polar'] / self.sequence_length * 100.0) + '\n' + \
               '    Charged\t\t' + str(self.amino_acids_composition['Charged']) + '\t\t\t' + str(
            self.amino_acids_composition['Charged'] / self.sequence_length * 100.0) + '\n' + \
               '    Basic\t\t' + str(self.amino_acids_composition['Basic']) + '\t\t\t' + str(
            self.amino_acids_composition['Basic'] / self.sequence_length * 100.0) + '\n' + \
               '    Acidic\t\t' + str(self.amino_acids_composition['Acidic']) + '\t\t\t' + str(
            self.amino_acids_composition['Acidic'] / self.sequence_length * 100.0) + '\n' + \
               '  Aromaticity: ' + str(self.aromaticity) + '\n' + \
               '  Instability: ' + str(self.instability) + '\n' + \
               '  Flexibility: ' + str(self.flexibility) + '\n' + \
               '  Weight list: ' + str(self.weight_list) + '\n' + \
               '  Protein scale: ' + str(self.protein_scale) + '\n' + \
               '  Isoelectric point: ' + str(self.isoelectric_point) + '\n' + \
               '  Secondary structure fraction: ' + str(self.secondary_structure_fraction) + '\n' + \
               '  Molecular weight: ' + str(self.molecular_weight) + '\n' + \
               '  Kyte plot: ' + str(self.kyte_plot) + '\n' + \
               '  M: ' + str(self.M) + '\n' + \
               '  Z: ' + str(self.Z) + '\n'

    def __repr__(self):
        return 'Protein computational parameters:\n' + ProteinParameters.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.sequence_length == other.sequence_length and \
               self.amino_acids_percents == other.amino_acids_percents and \
               self.amino_acids_composition == other.amino_acids_composition and \
               self.aromaticity == other.aromaticity and \
               self.instability == other.instability and \
               self.flexibility == other.flexibility and \
               self.weight_list == other.weight_list and \
               self.protein_scale == other.protein_scale and \
               self.isoelectric_point == other.isoelectric_point and \
               self.secondary_structure_fraction == other.secondary_structure_fraction and \
               self.molecular_weight == other.molecular_weight and \
               self.kyte_plot == other.kyte_plot and \
               self.M == other.M and \
               self.Z == other.Z


def count_acids_from_list(sequence, list):
    number = 0
    for acid in sequence:
        if acid in list:
            number += 1
    return number


def calculate_amino_acids_composition(sequence):
    composition = {}
    composition['Small'] = count_acids_from_list(sequence, ['A', 'B', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'])
    composition['Aliphatic'] = count_acids_from_list(sequence, ['A', 'I', 'L', 'V'])
    composition['Aromatic'] = count_acids_from_list(sequence, ['F', 'H', 'W', 'Y'])
    composition['Non-polar'] = count_acids_from_list(sequence, ['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'])
    composition['Polar'] = count_acids_from_list(sequence, ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z'])
    composition['Charged'] = count_acids_from_list(sequence, ['B', 'D', 'E', 'H', 'K', 'R', 'Z'])
    composition['Basic'] = count_acids_from_list(sequence, ['H', 'K', 'R'])
    composition['Acidic'] = count_acids_from_list(sequence, ['B', 'D', 'E', 'Z'])
    return composition
