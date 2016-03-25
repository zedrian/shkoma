from pandas import DataFrame, Series
from numpy import zeros, float64

from shkoma.utility import show_progress


def calculate_simple_statistics_for_serie(serie):
    described = serie.describe()

    stats = {}
    stats['mean'] = described['mean']
    stats['variance'] = serie.var()
    stats['skewness'] = serie.skew()
    stats['kurtosis'] = serie.kurtosis()

    stats['std'] = described['std']

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

    total_received_peptides_number = 0
    total_missed_peptides_number = 0
    for protein_record in protein_records:
        total_received_peptides_number += len(protein_record.received_peptide_records)
        total_missed_peptides_number += len(protein_record.missed_peptide_records)

    received_parameters = DataFrame(zeros((total_received_peptides_number, len(peptide_parameter_names)),
                                          dtype=float64), columns=peptide_parameter_names)
    missed_parameters = DataFrame(zeros((total_missed_peptides_number, len(peptide_parameter_names)),
                                        dtype=float64), columns=peptide_parameter_names)

    # fill received peptides parameters
    label = 'Filling received peptides parameter lists: '
    show_progress(label, 40, 0.0)
    index = 1
    for protein_record in protein_records:
        for received_peptide_record in protein_record.received_peptide_records:
            received_parameters['Sequence length'][index] = received_peptide_record.peptide_parameters.sequence_length
            received_parameters['Aromaticity'][index] = received_peptide_record.peptide_parameters.aromaticity
            received_parameters['Instability'][index] = received_peptide_record.peptide_parameters.instability
            received_parameters['Isoelectric point'][index] = \
                received_peptide_record.peptide_parameters.isoelectric_point
            received_parameters['Molecular weight'][index] = received_peptide_record.peptide_parameters.molecular_weight
            received_parameters['Kyte plot'][index] = received_peptide_record.peptide_parameters.kyte_plot
            received_parameters['Aliphatic index'][index] = received_peptide_record.peptide_parameters.aliphatic_index
            received_parameters['Boman index'][index] = received_peptide_record.peptide_parameters.boman_index
            received_parameters['Hydrophobicity'][index] = received_peptide_record.peptide_parameters.hydrophobicity

            for kidera_factor in received_peptide_record.peptide_parameters.kidera_factors:
                received_parameters['Kidera factor: {0}'.format(kidera_factor['name'])][index] = kidera_factor['value']

            show_progress(label, 40, index / total_received_peptides_number)
            index += 1
    print()

    # fill missed peptides parameters
    label = 'Filling missed peptides parameter lists: '
    show_progress(label, 40, 0.0)
    index = 1
    for protein_record in protein_records:
        for missed_peptide_record in protein_record.missed_peptide_records:
            missed_parameters['Sequence length'][index] = missed_peptide_record.peptide_parameters.sequence_length
            missed_parameters['Aromaticity'][index] = missed_peptide_record.peptide_parameters.aromaticity
            missed_parameters['Instability'][index] = missed_peptide_record.peptide_parameters.instability
            missed_parameters['Isoelectric point'][index] = missed_peptide_record.peptide_parameters.isoelectric_point
            missed_parameters['Molecular weight'][index] = missed_peptide_record.peptide_parameters.molecular_weight
            missed_parameters['Kyte plot'][index] = missed_peptide_record.peptide_parameters.kyte_plot
            missed_parameters['Aliphatic index'][index] = missed_peptide_record.peptide_parameters.aliphatic_index
            missed_parameters['Boman index'][index] = missed_peptide_record.peptide_parameters.boman_index
            missed_parameters['Hydrophobicity'][index] = missed_peptide_record.peptide_parameters.hydrophobicity

            for kidera_factor in missed_peptide_record.peptide_parameters.kidera_factors:
                missed_parameters['Kidera factor: {0}'.format(kidera_factor['name'])][index] = kidera_factor['value']

            show_progress(label, 40, index / total_missed_peptides_number)
            index += 1
    print()

    return received_parameters, missed_parameters


def fill_per_peptide_correlations(protein_records):
    kidera_factor_names = ['helix.bend.pref', 'side.chain.size', 'extended.str.pref',
                           'hydrophobicity', 'double.bend.pref', 'partial.spec.vol',
                           'flat.ext.pref', 'occurrence.alpha.reg', 'pK.C', 'surrounding.hydrop']

    per_peptide_correlation_parameter_names = ['Kidera factors', 'Amino acid percents', 'Amino acid compositions',
                                               'Charges', 'Hydrophobic moments']
    per_peptide_correlation_parameter_labels = ['{0} per peptide correlation (Pearson)'.format(name) for name in
                                                per_peptide_correlation_parameter_names]

    total_received_peptides_number = 0
    total_missed_peptides_number = 0
    for protein_record in protein_records:
        total_received_peptides_number += len(protein_record.received_peptide_records)
        total_missed_peptides_number += len(protein_record.missed_peptide_records)

    total_received_pairs_number = total_received_peptides_number * (total_received_peptides_number - 1) // 2
    received_per_peptide_correlations = DataFrame(zeros((total_received_pairs_number,
                                                         len(per_peptide_correlation_parameter_labels)),
                                                        dtype=float64),
                                                  columns=per_peptide_correlation_parameter_labels)
    total_missed_pairs_number = total_missed_peptides_number * (total_missed_peptides_number - 1) // 2
    missed_per_peptide_correlations = DataFrame(zeros((total_missed_pairs_number,
                                                       len(per_peptide_correlation_parameter_labels)),
                                                      dtype=float64),
                                                columns=per_peptide_correlation_parameter_labels)

    received_kidera_factors = DataFrame(zeros((len(kidera_factor_names), total_received_peptides_number),
                                              dtype=float64))
    missed_kidera_factors = DataFrame(zeros((len(kidera_factor_names), total_missed_peptides_number),
                                            dtype=float64))

    received_acid_percents = DataFrame(zeros((len('AGVMDYNSWLFIKPQCERTH'), total_received_peptides_number),
                                             dtype=float64))
    missed_acid_percents = DataFrame(zeros((len('AGVMDYNSWLFIKPQCERTH'), total_missed_peptides_number),
                                           dtype=float64))

    acid_group_names = ['Small', 'Aliphatic', 'Aromatic',
                        'Non-polar', 'Polar', 'Charged',
                        'Basic', 'Acidic']
    received_acid_compounds = DataFrame(zeros((len(acid_group_names), total_received_peptides_number),
                                              dtype=float64))
    missed_acid_compounds = DataFrame(zeros((len(acid_group_names), total_missed_peptides_number),
                                            dtype=float64))

    # received_charges = []
    # missed_charges = []
    #
    # received_hydrophobic_moments = []
    # missed_hydrophobic_moments = []

    label = 'Filling received peptides array-like parameter lists: '
    show_progress(label, 35, 0.0)
    index = 1
    for protein_record in protein_records:
        for received_peptide_record in protein_record.received_peptide_records:
            kidera_factor_index = 0
            for kidera_factor in received_peptide_record.peptide_parameters.kidera_factors:
                received_kidera_factors[index - 1][kidera_factor_index] = kidera_factor['value']
                kidera_factor_index += 1

            acid_index = 0
            for acid in 'AGVMDYNSWLFIKPQCERTH':
                received_acid_percents[index - 1][acid_index] = \
                    received_peptide_record.peptide_parameters.amino_acid_percents[acid]
                acid_index += 1

            group_index = 0
            for group in received_peptide_record.peptide_parameters.amino_acids_composition:
                received_acid_compounds[index - 1][group_index] = group['percent']
                group_index += 1

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
            show_progress(label, 35, index / total_received_peptides_number)
            index += 1
    print()

    label = 'Filling missed peptides array-like parameter lists: '
    show_progress(label, 35, 0.0)
    index = 1
    for protein_record in protein_records:
        for missed_peptide_record in protein_record.missed_peptide_records:
            kidera_factor_index = 0
            for kidera_factor in missed_peptide_record.peptide_parameters.kidera_factors:
                missed_kidera_factors[index - 1][kidera_factor_index] = kidera_factor['value']
                kidera_factor_index += 1

            acid_index = 0
            for acid in 'AGVMDYNSWLFIKPQCERTH':
                missed_acid_percents[index - 1][acid_index] = \
                    missed_peptide_record.peptide_parameters.amino_acid_percents[acid]
                acid_index += 1

            group_index = 0
            for group in missed_peptide_record.peptide_parameters.amino_acids_composition:
                missed_acid_compounds[index - 1][group_index] = group['percent']
                group_index += 1

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
                show_progress(label, 35, index / total_missed_peptides_number)
                index += 1
    print()

    print('Calculating Kidera factors per peptide Pearson correlation (received peptides): ', end='')
    received_per_peptide_correlations['Kidera factors per peptide correlation (Pearson)'] = \
        convert_correlation_matrix_to_serie(received_kidera_factors.corr(method='pearson'), 'Kidera factors')
    print('done')

    print('Calculating Kidera factors per peptide Pearson correlation (missed peptides): ', end='')
    missed_per_peptide_correlations['Kidera factors per peptide correlation (Pearson)'] = \
        convert_correlation_matrix_to_serie(missed_kidera_factors.corr(method='pearson'), 'Kidera factors')
    print('done')

    print('Calculating amino acid percents per peptide Pearson correlation (received peptides): ', end='')
    received_per_peptide_correlations['Amino acid percents per peptide correlation (Pearson)'] = \
        convert_correlation_matrix_to_serie(received_acid_percents.corr(method='pearson'), 'Amino acid percents')
    print('done')

    print('Calculating amino acid percents per peptide Pearson correlation (missed peptides): ', end='')
    missed_per_peptide_correlations['Amino acid percents per peptide correlation (Pearson)'] = \
        convert_correlation_matrix_to_serie(missed_acid_percents.corr(method='pearson'), 'Amino acid percents')
    print('done')

    print('Calculating amino acid compositions per peptide Pearson correlation (received peptides): ', end='')
    received_per_peptide_correlations['Amino acid compositions per peptide correlation (Pearson)'] = \
        convert_correlation_matrix_to_serie(received_acid_compounds.corr(method='pearson'), 'Amino acid compositions')
    print('done')

    print('Calculating amino acid compositions per peptide Pearson correlation (missed peptides): ', end='')
    missed_per_peptide_correlations['Amino acid compositions per peptide correlation (Pearson)'] = \
        convert_correlation_matrix_to_serie(missed_acid_compounds.corr(method='pearson'), 'Amino acid compositions')
    print('done')

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
    return received_per_peptide_correlations, missed_per_peptide_correlations


def convert_correlation_matrix_to_serie(matrix, name):
    size = len(matrix)
    total_pairs_number = size * (size - 1) // 2

    serie = Series(zeros(total_pairs_number, dtype=float64), name=name)

    pair_index = 0
    for i in range(0, size - 1):
        for j in range(i + 1, size - 1):
            serie[pair_index] = matrix[i][j]
            pair_index += 1

    return serie


def calculate_simple_statistics(parameters, per_peptide_correlations=None):
    label = 'Calculating simple statistics: '
    show_progress(label, 40, 0.0)
    stats = {}

    total_stats_length = len(parameters.columns)
    if per_peptide_correlations is not None:
        total_stats_length += len(per_peptide_correlations.columns)

    index = 1
    for parameter_name in parameters.columns:
        stats[parameter_name] = calculate_simple_statistics_for_serie(parameters[parameter_name])
        show_progress(label, 40, index / total_stats_length)
        index += 1
    if per_peptide_correlations is not None:
        for parameter_name in per_peptide_correlations.columns:
            stats[parameter_name] = calculate_simple_statistics_for_serie(per_peptide_correlations[parameter_name])
            show_progress(label, 40, index / total_stats_length)
            index += 1
    print()

    return stats


def save_simple_statistics_to_csv(stats, file_name):
    label = 'Saving simple statistics to \'{0}\': '.format(file_name)
    show_progress(label, 40, 0.0)

    index = 1
    with open(file_name, 'w') as file:
        file.write('name;mean;variance;skewness;kurtosis;std\n')
        for name in stats:
            file.write('{0};{1};{2};{3};{4};{5}\n'.format(name, stats[name]['mean'], stats[name]['variance'],
                                                          stats[name]['skewness'], stats[name]['kurtosis'],
                                                          stats[name]['std']))
            show_progress(label, 40, index / len(stats))
            index += 1
    print()
