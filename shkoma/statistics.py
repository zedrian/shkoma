from pandas import DataFrame
from numpy import zeros, float64

from shkoma.utility import show_progress


def calculate_simple_statistics(serie):
    described = serie.describe()

    stats = {}
    stats['mean'] = described['mean']
    stats['variance'] = serie.var()
    stats['skewness'] = serie.skew()
    stats['kurtosis'] = serie.kurtosis()

    stats['std'] = described['std']
    stats['variation'] = 0  # TODO: check how to calculate it using pandas

    return stats


def fill_parameter_lists(protein_records):
    kidera_factor_names = ['helix.bend.pref', 'side.chain.size', 'extended.str.pref',
                           'hydrophobicity', 'double.bend.pref', 'partial.spec.vol',
                           'flat.ext.pref', 'occurrence.alpha.reg', 'pK.C', 'surrounding.hydrop']
    kidera_factor_labels = ['Kidera factor: {0}'.format(name) for name in kidera_factor_names]

    peptide_parameter_names = ['Sequence length', 'Aromaticity', 'Instability',
                               'Isoelectric point', 'Molecular weight', 'Kyte plot',
                               'Aliphatic index', 'Boman index', 'Hydrophobicity']
    for name in kidera_factor_labels:
        peptide_parameter_names.append(name)
    # received['Kidera factors per peptide correlation (Kendall)'] = []
    # received['Amino acid percents per peptide correlation (Kendall)'] = []
    # received['Amino acid compositions per peptide correlation (Kendall)'] = []
    # received['Charges per peptide correlation (Kendall)'] = []
    # received['Hydrophobic moments per peptide correlation (Kendall)'] = []

    total_received_peptides_number = 0
    total_missed_peptides_number = 0
    for protein_record in protein_records:
        total_received_peptides_number += len(protein_record.received_peptide_records)
        total_missed_peptides_number += len(protein_record.missed_peptide_records)

    received = DataFrame(zeros((total_received_peptides_number, len(peptide_parameter_names)), dtype=float64),
                         columns=peptide_parameter_names)
    missed = DataFrame(zeros((total_missed_peptides_number, len(peptide_parameter_names)), dtype=float64),
                       columns=peptide_parameter_names)
    # received_kidera_factors = []
    # missed_kidera_factors = []
    #
    # received_acid_percents = []
    # missed_acid_percents = []
    #
    # received_acid_compounds = []
    # missed_acid_compounds = []
    #
    # received_charges = []
    # missed_charges = []
    #
    # received_hydrophobic_moments = []
    # missed_hydrophobic_moments = []

    # fill received peptides parameters
    label = 'Filling received peptides parameter lists: '
    show_progress(label, 40, 0.0)
    index = 1
    for protein_record in protein_records:
        for received_peptide_record in protein_record.received_peptide_records:
            received['Sequence length'][index] = received_peptide_record.peptide_parameters.sequence_length
            received['Aromaticity'][index] = received_peptide_record.peptide_parameters.aromaticity
            received['Instability'][index] = received_peptide_record.peptide_parameters.instability
            received['Isoelectric point'][index] = received_peptide_record.peptide_parameters.isoelectric_point
            received['Molecular weight'][index] = received_peptide_record.peptide_parameters.molecular_weight
            received['Kyte plot'][index] = received_peptide_record.peptide_parameters.kyte_plot
            received['Aliphatic index'][index] = received_peptide_record.peptide_parameters.aliphatic_index
            received['Boman index'][index] = received_peptide_record.peptide_parameters.boman_index
            received['Hydrophobicity'][index] = received_peptide_record.peptide_parameters.hydrophobicity

            for kidera_factor in received_peptide_record.peptide_parameters.kidera_factors:
                received['Kidera factor: {0}'.format(kidera_factor['name'])][index] = kidera_factor['value']

            show_progress(label, 40, index / total_received_peptides_number)
            index += 1

            #     kidera_factors.append(kidera_factor['value'])
            # received_kidera_factors.append(kidera_factors)

            # acid_percents = []
            # for acid in 'AGVMDYNSWLFIKPQCERTH':
            #     acid_percents.append(received_peptide_record.peptide_parameters.amino_acid_percents[acid])
            # received_acid_percents.append(acid_percents)
            #
            # acid_compound = []
            # for group in received_peptide_record.peptide_parameters.amino_acids_composition:
            #     acid_compound.append(group['percent'])
            # received_acid_compounds.append(acid_compound)
            #
            # charges = []
            # for charge in received_peptide_record.peptide_parameters.charges:
            #     charges.append(charge['charge'])
            # received_charges.append(charges)
            #
            # hydrophobic_moments = []
            # for moment in received_peptide_record.peptide_parameters.hydrophobic_moments:
            #     if moment['name'] != 'Polygly-polypro helix':
            #         hydrophobic_moments.append(moment['moment'])
            # received_hydrophobic_moments.append(hydrophobic_moments)
            # kidera_factors = []
    print()

    # fill missed peptides parameters
    label = 'Filling missed peptides parameter lists: '
    show_progress(label, 40, 0.0)
    index = 1
    for protein_record in protein_records:
        for missed_peptide_record in protein_record.missed_peptide_records:
            missed['Sequence length'][index] = missed_peptide_record.peptide_parameters.sequence_length
            missed['Aromaticity'][index] = missed_peptide_record.peptide_parameters.aromaticity
            missed['Instability'][index] = missed_peptide_record.peptide_parameters.instability
            missed['Isoelectric point'][index] = missed_peptide_record.peptide_parameters.isoelectric_point
            missed['Molecular weight'][index] = missed_peptide_record.peptide_parameters.molecular_weight
            missed['Kyte plot'][index] = missed_peptide_record.peptide_parameters.kyte_plot
            missed['Aliphatic index'][index] = missed_peptide_record.peptide_parameters.aliphatic_index
            missed['Boman index'][index] = missed_peptide_record.peptide_parameters.boman_index
            missed['Hydrophobicity'][index] = missed_peptide_record.peptide_parameters.hydrophobicity

            for kidera_factor in missed_peptide_record.peptide_parameters.kidera_factors:
                missed['Kidera factor: {0}'.format(kidera_factor['name'])][index] = kidera_factor['value']

            show_progress(label, 40, index / total_missed_peptides_number)
            index += 1

            #     kidera_factors.append(kidera_factor['value'])
            # missed_kidera_factors.append(kidera_factors)
            #
            # acid_percents = []
            # for acid in 'AGVMDYNSWLFIKPQCERTH':
            #     acid_percents.append(missed_peptide_record.peptide_parameters.amino_acid_percents[acid])
            # missed_acid_percents.append(acid_percents)
            #
            # acid_compound = []
            # for group in missed_peptide_record.peptide_parameters.amino_acids_composition:
            #     acid_compound.append(group['percent'])
            # missed_acid_compounds.append(acid_compound)
            #
            # charges = []
            # for charge in missed_peptide_record.peptide_parameters.charges:
            #     charges.append(charge['charge'])
            # missed_charges.append(charges)
            #
            # hydrophobic_moments = []
            # for moment in missed_peptide_record.peptide_parameters.hydrophobic_moments:
            #     if moment['name'] != 'Polygly-polypro helix':
            #         hydrophobic_moments.append(moment['moment'])
            # missed_hydrophobic_moments.append(hydrophobic_moments)
    print()
    #
    # label = 'Calculating Kidera factors Kendall correlation (received peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_kidera in range(0, len(received_kidera_factors)):
    #     for second_kidera in range(first_kidera + 1, len(received_kidera_factors)):
    #         received['Kidera factors per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(received_kidera_factors[first_kidera],
    #                                   received_kidera_factors[second_kidera]).correlation)
    #     show_progress(label, 40, index / len(received_kidera_factors))
    #     index += 1
    # print()
    #
    # label = 'Calculating Kidera factors Kendall correlation (missed peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_kidera in range(0, len(missed_kidera_factors)):
    #     for second_kidera in range(first_kidera + 1, len(received_kidera_factors)):
    #         missed['Kidera factors per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(missed_kidera_factors[first_kidera],
    #                                   missed_kidera_factors[second_kidera]).correlation)
    #     show_progress(label, 40, index / len(missed_kidera_factors))
    #     index += 1
    # print()
    #
    # label = 'Calculating amino acid percents Kendall correlation (received peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_percent in range(0, len(received_acid_percents)):
    #     for second_percent in range(first_percent + 1, len(received_acid_percents)):
    #         received['Amino acid percents per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(received_acid_percents[first_percent],
    #                                   received_acid_percents[second_percent]).correlation)
    #     show_progress(label, 40, index / len(received_acid_percents))
    #     index += 1
    # print()
    #
    # label = 'Calculating amino acid percents Kendall correlation (missed peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_percent in range(0, len(missed_acid_percents)):
    #     for second_percent in range(first_percent + 1, len(missed_acid_percents)):
    #         missed['Amino acid percents per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(missed_acid_percents[first_percent],
    #                                   missed_acid_percents[second_percent]).correlation)
    #     show_progress(label, 40, index / len(missed_acid_percents))
    #     index += 1
    # print()
    #
    # label = 'Calculating amino acid compositions Kendall correlation (received peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_compound in range(0, len(received_acid_compounds)):
    #     for second_compound in range(first_compound + 1, len(received_acid_compounds)):
    #         received['Amino acid compositions per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(received_acid_compounds[first_compound],
    #                                   received_acid_compounds[second_compound]).correlation)
    #     show_progress(label, 40, index / len(received_acid_compounds))
    #     index += 1
    # print()
    #
    # label = 'Calculating amino acid compositions Kendall correlation (missed peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_compound in range(0, len(missed_acid_compounds)):
    #     for second_compound in range(first_compound + 1, len(missed_acid_compounds)):
    #         missed['Amino acid compositions per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(missed_acid_compounds[first_compound],
    #                                   missed_acid_compounds[second_compound]).correlation)
    #     show_progress(label, 40, index / len(missed_acid_compounds))
    #     index += 1
    # print()
    #
    # label = 'Calculating charges Kendall correlation (received peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_charges in range(0, len(received_charges)):
    #     for second_charges in range(first_charges + 1, len(received_charges)):
    #         received['Charges per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(received_charges[first_charges], received_charges[second_charges]).correlation)
    #     show_progress(label, 40, index / len(received_charges))
    #     index += 1
    # print()
    #
    # label = 'Calculating charges Kendall correlation (missed peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_charges in range(0, len(missed_charges)):
    #     for second_charges in range(first_charges + 1, len(missed_charges)):
    #         missed['Charges per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(missed_charges[first_charges], missed_charges[second_charges]).correlation)
    #     show_progress(label, 40, index / len(missed_charges))
    #     index += 1
    # print()
    #
    # label = 'Calculating hydrophobic moments Kendall correlation (received peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_moments in range(0, len(received_hydrophobic_moments)):
    #     for second_moments in range(first_moments + 1, len(received_hydrophobic_moments)):
    #         received['Hydrophobic moments per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(received_hydrophobic_moments[first_moments],
    #                                   received_hydrophobic_moments[second_moments]).correlation)
    #     show_progress(label, 40, index / len(received_hydrophobic_moments))
    #     index += 1
    #
    # label = 'Calculating hydrophobic moments Kendall correlation (missed peptides): '
    # show_progress(label, 40, 0.0)
    # index = 1
    # for first_moments in range(0, len(missed_hydrophobic_moments)):
    #     for second_moments in range(first_moments + 1, len(missed_hydrophobic_moments)):
    #         missed['Hydrophobic moments per peptide correlation (Kendall)'].append(
    #             statistics.kendalltau(missed_hydrophobic_moments[first_moments],
    #                                   missed_hydrophobic_moments[second_moments]).correlation)
    #     show_progress(label, 40, index / len(missed_hydrophobic_moments))
    #     index += 1

    return received, missed


def calculate_simple_statistics_for_parameters_list(list):
    label = 'Calculating simple statistics: '
    show_progress(label, 40, 0.0)
    stats = {}

    kidera_factor_names = ['helix.bend.pref', 'side.chain.size', 'extended.str.pref',
                           'hydrophobicity', 'double.bend.pref', 'partial.spec.vol',
                           'flat.ext.pref', 'occurrence.alpha.reg', 'pK.C', 'surrounding.hydrop']
    parameter_names = ['Sequence length', 'Aromaticity', 'Instability',
                       'Isoelectric point', 'Molecular weight', 'Kyte plot',
                       'Aliphatic index', 'Boman index', 'Hydrophobicity']
    for name in kidera_factor_names:
        parameter_names.append('Kidera factor: {0}'.format(name))
    # parameter_names.append('Kidera factors per peptide correlation (Kendall)')
    # parameter_names.append('Amino acid percents per peptide correlation (Kendall)')
    # parameter_names.append('Amino acid compositions per peptide correlation (Kendall)')
    # parameter_names.append('Charges per peptide correlation (Kendall)')
    # parameter_names.append('Hydrophobic moments per peptide correlation (Kendall)')

    index = 1
    for parameter_name in parameter_names:
        stats[parameter_name] = calculate_simple_statistics(list[parameter_name])
        show_progress(label, 40, index / len(parameter_names))
        index += 1
    print()

    return stats


def save_simple_statistics_to_csv(stats, file_name):
    label = 'Saving simple statistics to \'{0}\': '.format(file_name)
    show_progress(label, 40, 0.0)

    kidera_factor_names = ['helix.bend.pref', 'side.chain.size', 'extended.str.pref',
                           'hydrophobicity', 'double.bend.pref', 'partial.spec.vol',
                           'flat.ext.pref', 'occurrence.alpha.reg', 'pK.C', 'surrounding.hydrop']
    parameter_names = ['Sequence length', 'Aromaticity', 'Instability',
                       'Isoelectric point', 'Molecular weight', 'Kyte plot',
                       'Aliphatic index', 'Boman index', 'Hydrophobicity']
    for name in kidera_factor_names:
        parameter_names.append('Kidera factor: {0}'.format(name))
    # parameter_names.append('Kidera factors per peptide correlation (Kendall)')
    # parameter_names.append('Amino acid percents per peptide correlation (Kendall)')
    # parameter_names.append('Amino acid compositions per peptide correlation (Kendall)')
    # parameter_names.append('Charges per peptide correlation (Kendall)')
    # parameter_names.append('Hydrophobic moments per peptide correlation (Kendall)')

    index = 1
    with open(file_name, 'w') as file:
        file.write('name;mean;variance;skewness;kurtosis;std;variation\n')
        for name in parameter_names:
            file.write('{0};{1};{2};{3};{4};{5};{6}\n'.format(name, stats[name]['mean'], stats[name]['variance'],
                                                              stats[name]['skewness'], stats[name]['kurtosis'],
                                                              stats[name]['std'], stats[name]['variation']))
            show_progress(label, 40, index / len(parameter_names))
            index += 1
    print()
