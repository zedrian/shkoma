import uniprot
from numpy import genfromtxt

from shkoma.peptide import Peptide
from shkoma.peptide_match import PeptideMatch
from shkoma.peptide_record import PeptideRecord, find_peptide_record_with_peptide
from shkoma.protein import Protein, find_protein_with_id
from shkoma.protein_record import ProteinRecord, find_protein_record_with_protein
from shkoma.utility import b2str


# load all main data as table from .csv file
def load_main_data_from_csv(file_name):
    main_data = genfromtxt(file_name, dtype=None, delimiter=';', names=True)
    return main_data


# construct list of proteins using main data
def construct_proteins(main_data):
    proteins = []

    # 1. fill list with unique proteins
    for line in main_data:
        # 1.1. construct protein from current line
        current_protein = Protein(id=b2str(line['accession_number']), name=b2str(line['entry_name']),
                                  mw=line['protein_mw'], pI=line['protein_pI'])

        # 1.2. add if not already exists in list
        if current_protein not in proteins:
            proteins.append(current_protein)

    return proteins


# load sequences from uniprot.org and fill such field in instances of class Protein
def fill_protein_sequences(proteins):
    for i in range(0, len(proteins)):
        print('Processing protein #' + str(i + 1) + ' of ' + str(len(proteins)) + '...')
        # 1. send requests while server will not get correct response
        server_response = None
        while not server_response:
            server_response = uniprot.get_metadata_with_some_seqid_conversions([proteins[i].id])

        # 2. store sequence in protein
        proteins[i].sequence = server_response[proteins[i].id]['sequence']


# save list of proteins to file
def save_proteins_to_csv(proteins, file_name):
    with open(file_name, 'w') as file:
        file.write('id;name;mw;pI;M;Z;sequence\n')
        for protein in proteins:
            file.write(protein.id + ';' + protein.name + ';' + str(protein.mw) + ';' + str(protein.pI) + ';' +
                       str(protein.M) + ';' + str(protein.Z) + ';' + str(protein.sequence) + '\n')


# load list of proteins from file
def load_proteins_from_csv(file_name):
    data = genfromtxt(file_name, dtype=None, delimiter=';', names=True)
    proteins = []
    for line in data:
        proteins.append(Protein(id=b2str(line['id']), name=b2str(line['name']), mw=line['mw'], pI=line['pI'], M=line['M'], Z=line['Z'], sequence=b2str(line['sequence'])))
    return proteins


# construct list of protein records using list of proteins and main data
def construct_protein_records(proteins, main_data):
    protein_records = []

    # 1. process all main data
    for line in main_data:
        # 1.1. construct peptide and peptide match from current analysis
        current_peptide = Peptide(sequence=b2str(line['sequence']), pI=line['peptide_pI'])
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

    return protein_records
