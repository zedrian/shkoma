class PeptideRecord:
    def __init__(self, peptide, matches=[]):
        self.peptide = peptide
        self.peptide_parameters = None
        self.matches = matches

    def __str__(self):
        if len(self.matches) != 0:
            return received_peptide_record_to_string(self)
        return missed_peptide_record_to_string(self)

    def __repr__(self):
        return '\nPeptide record:\n' + PeptideRecord.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.peptide == other.peptide and \
               self.peptide_parameters == other.peptide_parameters and \
               self.matches == other.matches


def find_peptide_record_with_peptide(records, peptide):
    for record in records:
        if record.peptide == peptide:
            return record

    return None


def received_peptide_record_to_string(record):
    # received peptide record interpretation = missed peptide record interpretation + matches
    result = missed_peptide_record_to_string(record)

    result += '  Matches: {0}\n'.format(len(record.matches))
    if len(record.matches) != 0:
        index = 1
        for match in record.matches:
            result += '  Match #{0} of {1}:\n'.format(index, len(record.matches)) + \
                      str(match)
            index += 1

    return result


def missed_peptide_record_to_string(record):
    result = str(record.peptide)
    if record.peptide_parameters is not None:
        result += str(record.peptide_parameters)
    return result
