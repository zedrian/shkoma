from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import *
from rpy2.robjects import r
from numpy import arange

from shkoma.additional_protein_parameter_data import *
from shkoma.enumerations import *


class PeptideParameters:
    def __init__(self, sequence):
        self.sequence = sequence
        self.sequence_length = len(sequence)
        analysis = ProteinAnalysis(sequence)

        self.amino_acid_percents = analysis.get_amino_acids_percent()
        self.amino_acids_composition = calculate_amino_acids_composition(sequence)
        self.aromaticity = analysis.aromaticity()
        self.instability = analysis.instability_index()
        self.flexibility = calculate_flexibility(sequence)
        protein_scale_parameters = [{'name': 'Hydrophilicity', 'dictionary': hw},
                                    {'name': 'Surface accessibility', 'dictionary': em},
                                    {'name': 'Janin Interior to surface transfer energy scale', 'dictionary': ja},
                                    {'name': 'Bulkiness', 'dictionary': bulkiness},
                                    {'name': 'Polarity', 'dictionary': polarity},
                                    {'name': 'Buried residues', 'dictionary': buried_residues},
                                    {'name': 'Average area buried', 'dictionary': average_area_buried},
                                    {'name': 'Retention time', 'dictionary': retention_time}]
        self.protein_scales = calculate_protein_scales(analysis, protein_scale_parameters)
        self.isoelectric_point = analysis.isoelectric_point()
        self.secondary_structure_fraction = calculate_secondary_structure_fraction(analysis)
        self.molecular_weight = analysis.molecular_weight()
        self.kyte_plot = analysis.gravy()
        self.pefing = calculate_pefing(sequence)

        # next parameters are calculated using R.Peptides
        r('require(Peptides)')
        r('sequence = "{0}"'.format(sequence))
        self.aliphatic_index = r('aindex(sequence)')[0]
        self.boman_index = r('boman(sequence)')[0]
        self.charges = calculate_charges(sequence, 1.0, 14.0, 0.5, 'Lehninger')
        self.hydrophobicity = r('seq(sequence)')[0]
        angles = [{'name': 'Alpha-helix', 'angle': -47},
                  {'name': '3-10-helix', 'angle': -26},
                  {'name': 'Pi-helix', 'angle': -80},
                  {'name': 'Omega', 'angle': 180},
                  {'name': 'Antiparallel beta-sheet', 'angle': 135},
                  {'name': 'Parallel beta-sheet', 'angle': 113}]
        if self.amino_acid_percents['P'] + self.amino_acid_percents['G'] > 0.3:
            angles.append({'name': 'Polygly-polypro helix', 'angle': 153})
        self.hydrophobic_moments = calculate_hydrophobic_moments(sequence, angles)
        self.kidera_factors = calculate_kidera_factors(sequence)
        self.peptide_types = calculate_peptide_types(sequence, angles)

    def __str__(self):
        return '  Sequence length: ' + str(self.sequence_length) + '\n' + \
               amino_acids_percents_to_string(self.amino_acid_percents, '  ') + \
               amino_acids_composition_to_string(self.amino_acids_composition, '  ') + \
               '  Aromaticity: {0:.3f}\n'.format(self.aromaticity) + \
               '  Instability: {0:.3f}\n'.format(self.instability) + \
               '  Flexibility: ' + str(self.flexibility) + '\n' + \
               protein_scales_to_string(self.protein_scales, '  ') + \
               '  Isoelectric point: {0:.3f}\n'.format(self.isoelectric_point) + \
               secondary_structure_fraction_to_string(self.secondary_structure_fraction, '  ') + \
               '  Molecular weight: {0:.3f}\n'.format(self.molecular_weight) + \
               '  Kyte plot: {0:.3f}\n'.format(self.kyte_plot) + \
               pefing_to_string(self.pefing, '  ') + \
               '  Aliphatic_index: {0:.3f}\n'.format(self.aliphatic_index) + \
               '  Boman index: {0:.3f}\n'.format(self.boman_index) + \
               charges_to_string(self.charges, '  ') + \
               '  Hydrophobicity: {0:.3f}\n'.format(self.hydrophobicity) + \
               hydrophobic_moments_to_string(self.hydrophobic_moments, '  ') + \
               kidera_factors_to_string(self.kidera_factors, '  ') + \
               peptide_types_to_string(self.peptide_types, '  ')

    def __repr__(self):
        return 'Peptide computational parameters:\n' + PeptideParameters.__str__(self)

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
    composition = []
    acid_groups = [{'name': 'Small', 'acids': ['A', 'B', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V']},
                   {'name': 'Aliphatic', 'acids': ['A', 'I', 'L', 'V']},
                   {'name': 'Aromatic', 'acids': ['F', 'H', 'W', 'Y']},
                   {'name': 'Non-polar', 'acids': ['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y']},
                   {'name': 'Polar', 'acids': ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z']},
                   {'name': 'Charged', 'acids': ['B', 'D', 'E', 'H', 'K', 'R', 'Z']},
                   {'name': 'Basic', 'acids': ['H', 'K', 'R']},
                   {'name': 'Acidic', 'acids': ['B', 'D', 'E', 'Z']}]

    for group in acid_groups:
        current_group_composition = {}
        current_group_composition['name'] = group['name']
        current_group_composition['amount'] = count_acids_from_list(sequence, group['acids'])
        current_group_composition['percent'] = current_group_composition['amount'] / len(sequence)
        composition.append(current_group_composition)

    return composition


def calculate_flexibility(sequence):
    if len(sequence) == 1:
        return [Flex[sequence[0]]]

    flexibilities = Flex
    window_size = 9
    weights = [0.25, 0.4375, 0.625, 0.8125, 1]
    sum = 5.25
    if len(sequence) < window_size:
        window_size = len(sequence)

        unit = 2 * 0.6 / (window_size - 1)
        weights = [0.0] * (window_size // 2)
        sum = 0
        for i in range(window_size // 2):
            weights[i] = 0.4 + unit * i
            sum += weights[i]
        sum *= 2
        if window_size % 2 == 1:
            sum += 1
            weights.append(1.0)

    scores = []

    for i in range(len(sequence) - window_size):
        subsequence = sequence[i:i + window_size]
        score = 0.0

        for j in range(window_size // 2):
            front = subsequence[j]
            back = subsequence[window_size - j - 1]
            score += (flexibilities[front] + flexibilities[back]) * weights[j]

        middle = subsequence[window_size // 2 + 1]
        score += flexibilities[middle]

        scores.append(score / sum)

    return scores


def calculate_protein_scales(analysis, protein_scale_parameters):
    protein_scales = []

    window_size = 9
    if analysis.length < 9:
        window_size = analysis.length
    for parameter in protein_scale_parameters:
        scale = {}
        scale['name'] = parameter['name']
        if window_size == 1:
            scale['value'] = None
        else:
            scale['value'] = analysis.protein_scale(parameter['dictionary'], window=window_size, edge=1.0)
        protein_scales.append(scale)

    return protein_scales


def protein_scales_to_string(protein_scales, prefix):
    result = prefix + 'Protein scales:\n'

    for scale in protein_scales:
        result += prefix + '  {0}: {1}\n'.format(scale['name'], scale['value'])

    return result


def calculate_secondary_structure_fraction(analysis):
    fraction = analysis.secondary_structure_fraction()

    return [{'name': 'Helix', 'value': fraction[0]},
            {'name': 'Turn', 'value': fraction[1]},
            {'name': 'Sheet', 'value': fraction[2]}]


def secondary_structure_fraction_to_string(fraction, prefix):
    result = prefix + 'Secondary structure fraction:\n'
    for shape in fraction:
        result += prefix + '  {0}: {1:.3f}\n'.format(shape['name'], shape['value'])
    return result


def amino_acids_percents_to_string(percents, prefix):
    result = prefix + 'Amino acid percents:\n' + \
             prefix + '  {0:<7}{1:>10}\n'.format('ACID', 'PERCENT')

    acids = 'AGVMDYNSWLFIKPQCERTH'
    for acid in acids:
        result += prefix + '  {0:<7}{1:>10.3%}\n'.format(acid, percents[acid])

    return result


def amino_acids_composition_to_string(composition, prefix):
    result = prefix + 'Amino acids composition:\n' + \
             prefix + '  {0:<10}{1:>10}{2:>12}\n'.format('GROUP', 'NUMBER', 'PERCENT')

    for group in composition:
        result += prefix + '  {0:<10}{1:>10}{2:>12.3f}\n'.format(group['name'], group['amount'], group['percent'])

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


def calculate_hydrophobic_moments(sequence, angles):
    r('require(Peptides)')
    r('sequence = "{0}"'.format(sequence))

    moments = []

    for angle in angles:
        moment = {}
        moment['name'] = angle['name']
        moment['angle'] = angle['angle']
        moment['moment'] = r('hmoment(sequence, {0}, 9)'.format(angle['angle']))[0]
        moments.append(moment)

    return moments


def hydrophobic_moments_to_string(moments, prefix):
    result = prefix + 'Hydrophobic moments:\n'

    for moment in moments:
        result += prefix + '  {0} (angle = {1}): {2:.3f}\n'.format(moment['name'], moment['angle'], moment['moment'])

    return result


def calculate_kidera_factors(sequence):
    r('require(Peptides)')
    r('sequence = "{0}"'.format(sequence))

    factors = []
    for name in kidera_factor_names:
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
