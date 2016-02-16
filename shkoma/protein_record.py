from shkoma.protein import Protein


class ProteinRecord:
    def __init__(self, protein=Protein, received_peptides=[], missed_peptides=[]):
        self.protein = protein
        self.received_peptides = received_peptides
        self.missed_peptides = missed_peptides

    def __str__(self):
        return 'protein: ' + str(self.protein) + ', received_peptides: ' + str(self.received_peptides) + \
               ', missed_peptides: ' + str(self.missed_peptides)

    def __repr__(self):
        return 'protein record (' + ProteinRecord.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.protein == other.protein and \
               self.received_peptides == other.received_peptides and \
               self.missed_peptides == other.missed_peptides


def find_protein_record_with_protein(records, protein):
    for record in records:
        if record.protein == protein:
            return record

    return None
