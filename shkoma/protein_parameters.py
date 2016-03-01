from Bio.SeqUtils.ProtParam import ProteinAnalysis


class ProteinParameters:
    def __init__(self, sequence):
        self.sequence_length = len(sequence)
        analysis = ProteinAnalysis(sequence)

        self.amino_acid_percents = analysis.get_amino_acids_percent()
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
        self.pefing = calculate_pefing(sequence)

    def __str__(self):
        return '  Sequence length: ' + str(self.sequence_length) + '\n' + \
               amino_acids_percents_to_string(self.amino_acid_percents, '  ') + \
               amino_acids_composition_to_string(self.amino_acids_composition, '  ', self.sequence_length) + \
               '  Aromaticity: {0:.3f}\n'.format(self.aromaticity) + \
               '  Instability: {0:.3f}\n'.format(self.instability) + \
               '  Flexibility: ' + str(self.flexibility) + '\n' + \
               '  Weight list: ' + str(self.weight_list) + '\n' + \
               '  Protein scale: ' + str(self.protein_scale) + '\n' + \
               '  Isoelectric point: {0:.3f}\n'.format(self.isoelectric_point) + \
               '  Secondary structure fraction: ' + str(self.secondary_structure_fraction) + '\n' + \
               '  Molecular weight: {0:.3f}\n'.format(self.molecular_weight) + \
               '  Kyte plot: {0:.3f}\n'.format(self.kyte_plot) + \
               '  M: {0:.3f}\n'.format(self.M) + \
               '  Z: {0:.3f}\n'.format(self.Z) + \
               pefing_to_string(self.pefing, '  ')

    def __repr__(self):
        return 'Protein computational parameters:\n' + ProteinParameters.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.sequence_length == other.sequence_length and \
               self.amino_acid_percents == other.amino_acid_percents and \
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
               self.Z == other.Z and \
               self.pefing == other.pefing


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


def amino_acids_percents_to_string(percents, prefix):
    result = prefix + 'Amino acid percents:\n' + \
             prefix + '  {0:<7}{1:>10}\n'.format('ACID', 'PERCENT')

    acids = 'AGVMDYNSWLFIKPQCERTH'
    for acid in acids:
        result += prefix + '  {0:<7}{1:>10.3%}\n'.format(acid, percents[acid])

    return result


def amino_acids_composition_to_string(composition, prefix, sequence_length):
    result = prefix + 'Amino acids composition:\n' + \
             prefix + '  {0:<10}{1:>10}{2:>12}\n'.format('GROUP', 'NUMBER', 'PERCENT')

    groups = ['Small', 'Aliphatic', 'Aromatic', 'Non-polar', 'Polar', 'Charged', 'Basic', 'Acidic']
    for group in groups:
        result += prefix + '  {0:<10}{1:>10}{2:>12.3%}\n'.format(group, composition[group],
                                                                 composition[group] / sequence_length)

    return result


# calculate peptide fingerprint
def calculate_pefing(sequence):
    pefing = {}

    acids = 'AGVMDYNSWLFIKPQCERTH'

    for first_acid in acids:
        line = {}
        first_acid_number = 0
        for second_acid in acids:
            line[second_acid] = 0
        for i in range(1, len(sequence)):
            if sequence[i-1] == first_acid:
                line[sequence[i]] += 1
                first_acid_number += 1

        if not first_acid_number == 0:
            for second_acid in acids:
                line[second_acid] /= first_acid_number
        pefing[first_acid] = line

    return pefing


def pefing_to_string(pefing, prefix):
    result = prefix + 'Peptide fingerprint:\n' + \
             prefix + '  {0:<5}'.format(' ')

    acids = 'AGVMDYNSWLFIKPQCERTH'

    for acid in acids:
        result += '{0:>10}'.format(acid)
    result += '\n'

    for first_acid in acids:
        result += prefix + '  {0:<5}'.format(first_acid)
        for second_acid in acids:
            result += '{0:>10.3%}'.format(pefing[first_acid][second_acid])
        result += '\n'

    return result
