from Bio.SeqUtils.ProtParam import ProteinAnalysis
from rpy2.robjects import r
from numpy import arange


class ProteinParameters:  # TODO: change to PeptideParameters
    def __init__(self, sequence):
        self.sequence = sequence
        self.sequence_length = len(sequence)
        analysis = ProteinAnalysis(sequence)

        self.amino_acid_percents = analysis.get_amino_acids_percent()
        self.amino_acids_composition = calculate_amino_acids_composition(sequence)
        self.aromaticity = analysis.aromaticity()
        self.instability = analysis.instability_index()
        self.flexibility = analysis.flexibility()
        self.weight_list = None  # analysis.weight_list(11, ?)  # TODO: remove
        self.protein_scale = None  # analysis.protein_scale(?, 11, ?)
        self.isoelectric_point = analysis.isoelectric_point()
        self.secondary_structure_fraction = analysis.secondary_structure_fraction()
        self.molecular_weight = analysis.molecular_weight()
        self.kyte_plot = analysis.gravy()
        self.M = 0  # TODO: move to ProteinParameters
        self.Z = 0  # TODO: move to ProteinParameters
        self.pefing = calculate_pefing(sequence)

        # next parameters are calculated using R.Peptides
        r('require(Peptides)')
        r('sequence = "{0}"'.format(sequence))
        self.aliphatic_index = r('aindex(sequence)')[0]
        self.boman_index = r('boman(sequence)')[0]
        self.charges = calculate_charges(sequence, 1.0, 14.0, 0.5, 'Lehninger')
        self.hydrophobicity = r('seq(sequence)')[0]
        self.hydrophobic_moments = calculate_hydrophobic_moments(sequence)
        self.kidera_factors = calculate_kidera_factors(sequence)
        angles = [{'name': 'Alpha-helix', 'angle': -47},
                  {'name': '3-10-helix', 'angle': -26},
                  {'name': 'Pi-helix', 'angle': -80},
                  {'name': 'Omega', 'angle': 180},
                  {'name': 'Antiparallel beta-sheet', 'angle': 135},
                  {'name': 'Parallel beta-sheet', 'angle': 113}]
        self.peptide_types = calculate_peptide_types(sequence, angles)

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
               pefing_to_string(self.pefing, '  ') + \
               '  Aliphatic_index: {0:.3f}\n'.format(self.aliphatic_index) + \
               '  Boman index: {0:.3f}\n'.format(self.boman_index) + \
               charges_to_string(self.charges, '  ') + \
               '  Hydrophobicity: {0:.3f}\n'.format(self.hydrophobicity) + \
               hydrophobic_moments_to_string(self.hydrophobic_moments, '  ') + \
               kidera_factors_to_string(self.kidera_factors, '  ') + \
               peptide_types_to_string(self.peptide_types, '  ')

    def __repr__(self):
        return 'Protein computational parameters:\n' + ProteinParameters.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.sequence == other.sequence


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
            if sequence[i - 1] == first_acid:
                line[sequence[i]] += 1
                first_acid_number += 1

        if not first_acid_number == 0:
            for second_acid in acids:
                line[second_acid] /= first_acid_number
        pefing[first_acid] = line

    return pefing


def pefing_to_string(pefing, prefix):
    result = prefix + 'Peptide fingerprint:\n' + \
             prefix + '  {0:<2}'.format(' ')

    acids = 'AGVMDYNSWLFIKPQCERTH'

    for acid in acids:
        result += '{0:^10}'.format(acid)
    result += '\n'

    for first_acid in acids:
        result += prefix + '  {0:<2}'.format(first_acid)
        for second_acid in acids:
            if pefing[first_acid][second_acid] > 0:
                result += '{0:>10.3%}'.format(pefing[first_acid][second_acid])
            else:
                result += '{0:>10}'.format('.')
        result += '\n'

    return result


def calculate_charges(sequence, ph_min, ph_max, ph_step, method):
    r('require(Peptides)')
    r('sequence = "{0}"'.format(sequence))

    charges = []
    for ph in arange(ph_min, ph_max + ph_step, ph_step):
        charge = {}
        charge['pH'] = ph
        charge['charge'] = r('charge(sequence, {0}, "{1}")'.format(ph, method))[0]
        charges.append(charge)

    return charges


def charges_to_string(charges, prefix):
    result = prefix + 'Charges for different pH:\n' + \
             prefix + '  {0:>4}{1:>10}\n'.format('pH', 'CHARGE')

    for charge in charges:
        result += prefix + '  {0:>4.1f}{1:>10.3f}\n'.format(charge['pH'], charge['charge'])

    return result


def calculate_hydrophobic_moments(sequence):
    r('require(Peptides)')
    r('sequence = "{0}"'.format(sequence))

    moment100 = {}
    moment100['shape'] = 'Helix'
    moment100['moment'] = r('hmoment(sequence, 100, 11)')[0]

    moment160 = {}
    moment160['shape'] = 'Beta-sheet'
    moment160['moment'] = r('hmoment(sequence, 160, 11)')[0]

    return [moment100, moment160]


def hydrophobic_moments_to_string(moments, prefix):
    result = prefix + 'Hydrophobic moments:\n'

    for moment in moments:
        result += prefix + '  {0}: {1:.3f}\n'.format(moment['shape'], moment['moment'])

    return result


def calculate_kidera_factors(sequence):
    r('require(Peptides)')
    r('sequence = "{0}"'.format(sequence))

    names = ['helix.bend.pref', 'side.chain.size', 'extended.str.pref',
             'hydrophobicity', 'double.bend.pref', 'partial.spec.vol',
             'flat.ext.pref', 'occurrence.alpha.reg', 'pK.C', 'surrounding.hydrop']

    factors = []
    for name in names:
        factor = {}
        factor['name'] = name
        factor['value'] = r('kidera(seq=sequence, factor="{0}")'.format(name))[0]
        factors.append(factor)

    return factors


def kidera_factors_to_string(factors, prefix):
    result = prefix + 'Kidera factors:\n'

    for factor in factors:
        result += prefix + '  {0}: {1:.3f}\n'.format(factor['name'], factor['value'])

    return result


def calculate_peptide_types(sequence, angles):
    r('require(Peptides)')
    r('sequence = "{0}"'.format(sequence))

    types = []
    for angle in angles:
        data = r('membpos(sequence, {0})'.format(angle['angle']))
        current_types = {}
        current_types['name'] = angle['name']
        current_types['angle'] = angle['angle']
        current_types['types'] = []
        for i in range(0, len(data[0])):
            line = {}
            line['peptide'] = data[0][i]
            line['H'] = data[1][i]
            line['uH'] = data[2][i]
            line['memb position'] = data[3][i]
            current_types['types'].append(line)
        types.append(current_types)

    return types


def peptide_types_to_string(types, prefix):
    result = prefix + 'Peptide types, calculated for different angles:\n'

    for current_types in types:
        result += prefix + '  {0} (angle = {1}):\n'.format(current_types['name'], current_types['angle']) + \
                  prefix + '    {0:>15}{1:>10}{2:>10}{3:>15}\n'.format('PEPTIDE', 'H', 'uH', 'MEMB_POSITION')
        for line in current_types['types']:
            result += prefix + '    {0:>15}{1:>10.3f}{2:>10.3f}{3:>15}\n'.format(line['peptide'], line['H'], line['uH'],
                                                                                 line['memb position'])

    return result
