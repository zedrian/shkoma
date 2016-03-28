from shkoma import correlation

main_data = correlation.load_main_data_from_csv(raw_input('Please, enter the .scv file name:' + '\n'))

proteins = correlation.construct_proteins(main_data)
correlation.fill_protein_sequences(proteins)

user_choice = raw_input('Do you want to save your proteins into a new file? "Yes/No"')
if user_choice.toUpper() == 'YES':
    protein_file_name = raw_input('Please, enter new file name')
    correlation.save_proteins_to_csv(proteins, protein_file_name + '.csv')

protein_records = correlation.construct_protein_records(proteins, main_data)

correlation.fill_missed_peptide_records(protein_records)

correlation.fill_peptide_parameters(protein_records)

user_choice = raw_input('Do you want to save your proteins into a new file? "Yes/No"')
if user_choice.toUpper() == 'Yes':
    folder_name = raw_input('Please, enter new folder name')
    correlation.save_protein_records_to_folder(protein_records, folder_name)