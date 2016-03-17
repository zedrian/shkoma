if __name__ == '__main__':
    from shkoma import correlation

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

    protein_records = correlation.construct_protein_records(proteins, main_data)
    correlation.fill_missed_peptide_records(protein_records)
    correlation.fill_peptide_parameters(protein_records)

    user_choice = input('Do you want to save your protein records [Yes/No]? ')
    if user_choice.upper() == 'YES':
        folder_name = input('Please, enter folder name: ')
        correlation.save_protein_records_to_folder(protein_records, folder_name)
