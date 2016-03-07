from shkoma.peptide_record import received_peptide_record_to_string, missed_peptide_record_to_string


class ProteinRecord:
    def __init__(self, protein, protein_parameters=None, received_peptide_records=[], missed_peptide_records=[]):
        self.protein = protein
        self.protein_parameters = protein_parameters
        self.received_peptide_records = received_peptide_records
        self.missed_peptide_records = missed_peptide_records

    def __str__(self):
        result = str(self.protein)
        if self.protein_parameters is not None:
            result += str(self.protein_parameters)

        result += '\nReceived peptide records: {0}\n'.format(len(self.received_peptide_records))
        index = 1
        for record in self.received_peptide_records:
            result += 'Peptide record #{0} of {1}:\n'.format(index, len(self.received_peptide_records)) + \
                      '{0}\n'.format(received_peptide_record_to_string(record))
            index += 1

        result += '\nMissed peptide records: {0}\n'.format(len(self.missed_peptide_records))
        index = 1
        for record in self.missed_peptide_records:
            result += 'Peptide record #{0} of {1}:\n'.format(index, len(self.missed_peptide_records)) + \
                      '{0}\n'.format(missed_peptide_record_to_string(record))
            index += 1

        return result

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
