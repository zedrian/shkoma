from shkoma.protein import Protein


class ProteinRecord:
    def __init__(self, protein, received_peptide_records=[], missed_peptide_records=[]):
        self.protein = protein
        self.received_peptide_records = received_peptide_records
        self.missed_peptide_records = missed_peptide_records

    def __str__(self):
        return 'Protein record: ' + str(self.protein) + \
               ', received_peptide_records: ' + str(self.received_peptide_records) + \
               ', missed_peptide_records: ' + str(self.missed_peptide_records)

    def __repr__(self):
        return 'protein record (' + ProteinRecord.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.protein == other.protein and \
               self.received_peptide_records == other.received_peptide_records and \
               self.missed_peptide_records == other.missed_peptide_records


def find_protein_record_with_protein(records, protein):
    for record in records:
        if record.protein == protein:
            return record

    return None
