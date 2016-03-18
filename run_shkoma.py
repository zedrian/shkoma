if __name__ == '__main__':
    from shkoma import correlation, statistics

    main_data = correlation.load_main_data_from_csv(input('Please, enter the main data .scv file name: '))

    user_choice = input('Do you want to construct proteins using main data [Yes] or to load them from file [No]? ')
    if user_choice.upper() == 'YES':
        proteins = correlation.construct_proteins(main_data)
        correlation.fill_protein_sequences(proteins)

        user_choice = input('Do you want to save your proteins into a new file [Yes/No]? ')
        if user_choice.upper() == 'YES':
            protein_file_name = input('Please, enter new file name: ')
            correlation.save_proteins_to_csv(proteins, protein_file_name + '.csv')
    else:
        proteins = correlation.load_proteins_from_csv(input('Please, enter proteins file name: '))

    protein_records = correlation.construct_protein_records(proteins, main_data)[1:3]
    correlation.fill_missed_peptide_records(protein_records)
    correlation.fill_peptide_parameters(protein_records)

    user_choice = input('Do you want to save your protein records [Yes/No]? ')
    if user_choice.upper() == 'YES':
        folder_name = input('Please, enter folder name: ')
        correlation.save_protein_records_to_folder(protein_records, folder_name)

    received_parameters, missed_parameters = statistics.fill_parameter_lists(protein_records)
    received_per_peptide_correlations, missed_per_peptide_correlations = statistics.fill_per_peptide_correlations(protein_records)
    received_statistics = statistics.calculate_simple_statistics(received_parameters, received_per_peptide_correlations)
    missed_statistics = statistics.calculate_simple_statistics(missed_parameters, missed_per_peptide_correlations)

    user_choice = input('Do you want to save statistics for received peptides [Yes/No]? ')
    if user_choice.upper() == 'YES':
        file_name = input('Please, enter new file name: ')
        statistics.save_simple_statistics_to_csv(received_statistics, file_name)

    user_choice = input('Do you want to save statistics for missed peptides [Yes/No]? ')
    if user_choice.upper() == 'YES':
        file_name = input('Please, enter new file name: ')
        statistics.save_simple_statistics_to_csv(missed_statistics, file_name)
