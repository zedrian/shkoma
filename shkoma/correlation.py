import uniprot
from numpy import genfromtxt

from shkoma.protein import Protein


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
        print('processing protein #' + str(i + 1) + ' from ' + str(len(proteins)) + '...')
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