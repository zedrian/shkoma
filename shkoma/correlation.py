import uniprot
from numpy import genfromtxt

from shkoma.peptide import Peptide
from shkoma.protein import Protein, find_protein_with_id
from shkoma.protein_record import ProteinRecord, find_protein_record_with_protein


# load all main data as table from .csv file
def load_main_data_from_csv(file_name):
    main_data = genfromtxt(file_name, dtype=None, delimiter=';', names=True)
    return main_data


# simple function converting numpy._bytes to human-readable string
def b2str(bytes):
    return str(bytes)[2:-1]


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


# construct list of protein records using list of proteins and main data
def construct_protein_records(proteins, main_data):
    records = []

    # 1. process all main data
    for line in main_data:
        # 1.1. construct peptide from current analysis
        current_peptide = Peptide(sequence='some sequence')

        # 1.1. get protein id for current analysis
        current_protein_id = line['accession_number']

        # 1.2. find protein with such id
        protein = find_protein_with_id(proteins,
                                       current_protein_id)  # TODO: if such protein not exists, extract Protein object and add to proteins

        # 1.3. find record with such protein
        record = find_protein_record_with_protein(records, protein)

        # 1.4. if record with such protein exists, add current analysis as received_peptide
        if record is not None:
            record.received_peptides.append(current_peptide)
        # ...or create new protein record otherwise
        else:
            record = ProteinRecord(protein)
            record.received_peptides.append(current_peptide)
            records.append(record)

    return records
