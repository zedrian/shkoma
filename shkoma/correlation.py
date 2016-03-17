import uniprot
from numpy import genfromtxt
from sys import stdout

from shkoma.alignment import cut_received_peptide_sequences, trypsinolize_sequence
from shkoma.peptide import Peptide
from shkoma.peptide_match import PeptideMatch
from shkoma.peptide_parameters import PeptideParameters
from shkoma.peptide_record import PeptideRecord, find_peptide_record_with_peptide
from shkoma.protein import Protein, find_protein_with_id
from shkoma.protein_record import ProteinRecord, find_protein_record_with_protein
from shkoma.protein_parameters import ProteinParameters
from shkoma.utility import b2str, show_progress


# load all main data as table from .csv file
def load_main_data_from_csv(file_name):
    label = 'Loading main data from \'{0}\': '.format(file_name)
    show_progress(label, 40, 0.0)

    main_data = genfromtxt(file_name, dtype=None, delimiter=';', names=True)
    show_progress(label, 40, 1.0)
    print()

    return main_data


# construct list of proteins using main data
def construct_proteins(main_data):
    proteins = []
    label = 'Constructing proteins from main data: '
    show_progress(label, 40, 0.0)

    # 1. fill list with unique proteins
    index = 1
    for line in main_data:
        # 1.1. construct protein from current line
        current_protein = Protein(id=b2str(line['accession_number']), name=b2str(line['entry_name']))

        # 1.2. add if not already exists in list
        if current_protein not in proteins:
            proteins.append(current_protein)

        show_progress(label, 40, index / len(main_data))
        index += 1
    print()

    return proteins


# load sequences from uniprot.org and fill such field in instances of class Protein
def fill_protein_sequences(proteins):
    label = 'Filling protein sequences: '

    show_progress(label, 40, 0.0)
    for i in range(0, len(proteins)):
        # 1. send requests while server will not get correct response
        server_response = None
        while not server_response:
            server_response = uniprot.get_metadata_with_some_seqid_conversions([proteins[i].id])

        # 2. store sequence in protein
        proteins[i].sequence = server_response[proteins[i].id]['sequence']

        show_progress(label, 40, i / len(proteins))
    print()


# save list of proteins to file
def save_proteins_to_csv(proteins, file_name):
    label = 'Saving proteins to \'{0}\': '.format(file_name)
    show_progress(label, 40, 0.0)

    with open(file_name, 'w') as file:
        file.write('id;name;sequence\n')
        index = 1
        for protein in proteins:
            file.write(protein.id + ';' + protein.name + ';' + protein.sequence + '\n')
            show_progress(label, 40, index / len(proteins))
            index += 1
    print()


# load list of proteins from file
def load_proteins_from_csv(file_name):
    label = 'Loading proteins from \'{0}\': '.format(file_name)
    show_progress(label, 40, 0.0)

    data = genfromtxt(file_name, dtype=None, delimiter=';', names=True)
    proteins = []
    index = 1
    for line in data:
        proteins.append(Protein(id=b2str(line['id']), name=b2str(line['name']), sequence=b2str(line['sequence'])))
        show_progress(label, 40, index / len(data))
        index += 1
    print()

    return proteins


# construct list of protein records and fill received peptide records using list of proteins and main data
def construct_protein_records(proteins, main_data):
    label = 'Constructing protein records: '
    show_progress(label, 40, 0.0)
    protein_records = []

    # 1. process all main data
    index = 1
    for line in main_data:
        # 1.1. construct peptide and peptide match from current analysis
        current_peptide = Peptide(sequence=b2str(line['sequence']))
        current_peptide_match = PeptideMatch(analysis_name=b2str(line['filename']), score=line['score'],
                                             reverse_score=line['reverseScore'],
                                             percent_of_scored_peak_intensity=line['percent_scored_peak_intensity'],
                                             total_intensity=line['totalIntensity'],
                                             precursor_averagine_chi_squared=line['precursorAveragineChiSquared'],
                                             retention_time_min=line['retentionTimeMin'],
                                             chromatographic_peak_width_in_seconds=line['chromatographicPeakWidthSec'])

        # 1.2. get protein id for current analysis
        current_protein_id = b2str(line['accession_number'])

        # 1.3. find protein with such id
        protein = find_protein_with_id(proteins,
                                       current_protein_id)  # TODO: if such protein not exists, extract Protein object and add to proteins

        # 1.4. find record with such protein
        protein_record = find_protein_record_with_protein(protein_records, protein)

        # 1.5. if record with such protein exists, add current match to received peptides
        if protein_record is not None:
            # 1.5.1. if such peptide was already received, add peptide match
            peptide_record = find_peptide_record_with_peptide(protein_record.received_peptide_records, current_peptide)
            if peptide_record is not None:
                peptide_record.matches.append(current_peptide_match)
            # 1.5.2. if such peptide was not received yet, add peptide record with this one peptide match
            else:
                current_peptide_record = PeptideRecord(current_peptide, [current_peptide_match])
                protein_record.received_peptide_records.append(current_peptide_record)
        # 1.6. if protein record with such protein not exists, create new protein record
        else:
            current_peptide_record = PeptideRecord(current_peptide, [current_peptide_match])
            protein_record = ProteinRecord(protein, received_peptide_records=[current_peptide_record])
            protein_records.append(protein_record)

        show_progress(label, 40, index / len(main_data))
        index += 1
    print()

    # 2. sort peptide records by length (starting from longest)
    label = 'Filling received peptide records: '
    show_progress(label, 40, 0.0)
    index = 1
    for protein_record in protein_records:
        protein_record.received_peptide_records = sorted(protein_record.received_peptide_records, key=lambda peptide_record: len(peptide_record.peptide.sequence), reverse=True)
        show_progress(label, 40, index / len(protein_records))
        index += 1
    print()

    return protein_records


# fill computational protein parameters for each protein record
def fill_protein_parameters(protein_records):
    label = 'Filling protein parameters: '
    show_progress(label, 40, 0.0)
    index = 1
    for protein_record in protein_records:
        protein_record.protein_parameters = ProteinParameters(protein_record.protein.sequence)
        show_progress(label, 40, index / len(protein_records))
        index += 1
    print()


# fill computational peptide parameters for each protein record
def fill_peptide_parameters(protein_records):
    protein_index = 1
    for protein_record in protein_records:
        print('Processing protein record #{0} of {1}:'.format(protein_index, len(protein_records)))
        stdout.flush()

        # 1. process received peptide records first
        label = '{0:>25}: '.format('Received peptides ({0})'.format(len(protein_record.received_peptide_records)))
        show_progress(label, 40, 0.0)
        peptide_index = 1
        for peptide_record in protein_record.received_peptide_records:
            peptide_record.peptide_parameters = PeptideParameters(peptide_record.peptide.sequence)
            show_progress(label, 40, peptide_index / len(protein_record.received_peptide_records))
            peptide_index += 1
        print()

        # 2. process then missed peptide records
        label = '{0:>25}: '.format('Missed peptides ({0})'.format(len(protein_record.missed_peptide_records)))
        show_progress(label, 40, 0.0)
        peptide_index = 1
        for peptide_record in protein_record.missed_peptide_records:
            peptide_record.peptide_parameters = PeptideParameters(peptide_record.peptide.sequence)
            show_progress(label, 40, peptide_index / len(protein_record.missed_peptide_records))
            peptide_index += 1
        print()

        protein_index += 1
        print()
    print('Processing protein records: done.')


# fill lists of missed peptide records for each protein record
def fill_missed_peptide_records(protein_records):
    label = 'Filling missed peptide records: '
    show_progress(label, 40, 0.0)
    index = 1
    for protein_record in protein_records:
        # 1. construct list of sequences of received peptides
        received_sequences = [peptide_record.peptide.sequence for peptide_record in protein_record.received_peptide_records]

        # 2. calculate list of missed sequence fragments
        missed_sequences = cut_received_peptide_sequences(protein_record.protein.sequence, received_sequences)
        missed_sequences = [trypsinolize_sequence(x) for x in missed_sequences]

        # 3. construct peptide record for each fragment and store them in missed peptide records
        protein_record.missed_peptide_records = []
        for missed_sequences_list in missed_sequences:
            for fragment in missed_sequences_list:
                protein_record.missed_peptide_records.append(PeptideRecord(peptide=Peptide(sequence=fragment)))

        show_progress(label, 40, index / len(protein_records))
        index += 1
    print()


# save list of protein records to folder (one file for each protein record)
def save_protein_records_to_folder(protein_records, folder='results/'):
    if not folder[-1] == '/':
        folder += '/'

    label = 'Saving protein records to \'{0}\': '.format(folder)
    show_progress(label, 40, 0.0)

    index = 1
    for protein_record in protein_records:
        with open(folder + protein_record.protein.id + '.txt', 'w') as file:
            file.write(str(protein_record))
        show_progress(label, 40, index / len(protein_records))
        index += 1
    print()
