class ProteinRecord:
    def __init__(self, protein, protein_parameters=None, received_peptide_records=[], missed_peptide_records=[]):
        self.protein = protein
        self.protein_parameters = protein_parameters
        self.received_peptide_records = received_peptide_records
        self.missed_peptide_records = missed_peptide_records

    def __str__(self):
        return str(self.protein) + \
               str(self.protein_parameters) + \
               'Received peptide records: ' + str(len(self.received_peptide_records)) + '\n' + \
               str(self.received_peptide_records) + '\n' + \
               'Missed peptide records: ' + str(len(self.missed_peptide_records)) + '\n' + \
               str(self.missed_peptide_records) + '\n'

    def __repr__(self):
        return '\nProtein record:\n' + ProteinRecord.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.protein == other.protein and \
               self.protein_parameters == other.protein_parameters and \
               self.received_peptide_records == other.received_peptide_records and \
               self.missed_peptide_records == other.missed_peptide_records


def find_protein_record_with_protein(records, protein):
    for record in records:
        if record.protein == protein:
            return record

    return None
