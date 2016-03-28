from shkoma import correlation

protein_file_boolean = False
protein_parameters_boolean = False

main_data = correlation.load_main_data_from_csv(input('Please, enter the .scv file name:' + '\n'))
user_choice_protein = input('Do you want to save your proteins into a new file? "Y/N"' + '\n')
if user_choice_protein.toUpper() == 'Y':
    protein_file_name = input('Please, enter new file name' + '\n')
    protein_file_boolean = True
user_choice_folder = input('Do you want to save your protein parameters into a new folder? "Y/N"' + '\n')
if user_choice_folder.toUpper() == 'Y':
    folder_name = input('Please, enter new folder name' + '\n')
    protein_parameters_boolean = True

proteins = correlation.construct_proteins(main_data)
correlation.fill_protein_sequences(proteins)

if protein_file_boolean:
    correlation.save_proteins_to_csv(proteins, protein_file_name + '.csv')

protein_records = correlation.construct_protein_records(proteins, main_data)

correlation.fill_missed_peptide_records(protein_records)

correlation.fill_peptide_parameters(protein_records)

if protein_parameters_boolean:
    correlation.save_protein_records_to_folder(protein_records, folder_name)