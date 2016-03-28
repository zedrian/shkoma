kidera_factor_names = ['helix.bend.pref', 'side.chain.size', 'extended.str.pref',
                       'hydrophobicity', 'double.bend.pref', 'partial.spec.vol',
                       'flat.ext.pref', 'occurrence.alpha.reg', 'pK.C', 'surrounding.hydrop']
kidera_factor_labels = ['Kidera factor: {0}'.format(name) for name in kidera_factor_names]


amino_acid_group_names = ['Small', 'Aliphatic', 'Aromatic',
                          'Non-polar', 'Polar', 'Charged',
                          'Basic', 'Acidic']

hydrophobic_moments_names = ['Alpha-helix', '3-10-helix', 'Pi-helix',
                             'Omega', 'Antiparallel beta-sheet', 'Parallel beta-sheet']

peptide_parameter_names = ['Sequence length', 'Aromaticity', 'Instability',
                           'Isoelectric point', 'Molecular weight', 'Kyte plot',
                           'Aliphatic index', 'Boman index', 'Hydrophobicity']
for name in kidera_factor_labels:
    peptide_parameter_names.append(name)

per_peptide_correlation_parameter_names = ['Kidera factors', 'Amino acid percents', 'Amino acid compositions',
                                           'Charges', 'Hydrophobic moments', 'Secondary structure fractions']