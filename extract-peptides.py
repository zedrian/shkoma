from os import listdir
from os.path import isfile, join
from sys import argv

from shkoma import correlation, statistics


def print_usage():
    print('----------------------------------------------------------------------')
    print('command line usage: ')
    print('  extract-peptides.py [experiments_folder] [peptides_folder]')
    print('where')
    print('  experiments_folder - folder with .csv`s with your experimental data,')
    print('  peptides_folder - folder where to put .csv`s with peptide parameters')
    print('')
    print('output:')
    print('  [experiment].peptides.csv and [experiment].statistics.csv,')
    print('  where [experiment] - file name from experiments_folder')
    print('----------------------------------------------------------------------')


def main():
    if len(argv) == 1:
        print_usage()
        experiments_folder = input('Enter folder name with experimental data: ')
        results_folder = input('Enter folder name for results: ')
    else:
        experiments_folder = argv[1]
        results_folder = argv[2]
    experiments_files = [f for f in listdir(experiments_folder) if isfile(join( experiments_folder, f))]

    for experiment_file in experiments_files:
        experiment_file_name = join(experiments_folder, experiment_file)
        main_data = correlation.load_main_data_from_csv(experiment_file_name)
        proteins = correlation.construct_proteins(main_data)
        protein_records = correlation.construct_protein_records(proteins, main_data)
        correlation.fill_peptide_parameters(protein_records)

        received, missed = statistics.fill_parameter_lists(protein_records)
        peptides_file_name = join(results_folder, experiment_file[:-4] + '.peptides.csv')
        statistics.save_peptide_parameter_lists_to_csv(received, peptides_file_name)

        stats = statistics.calculate_simple_statistics(received)
        statistics_file_name = join(results_folder, experiment_file[:-4] + '.statistics.csv')
        statistics.save_simple_statistics_to_csv(stats, statistics_file_name)


if __name__ == '__main__':
    main()
