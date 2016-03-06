class Peptide:
    def __init__(self, sequence):
        self.sequence = sequence

    def __str__(self):
        return '  Sequence: {0}\n'.format(self.sequence)

    def __repr__(self):
        return 'Peptide general information:\n' + Peptide.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.sequence == other.sequence  # TODO: add other checks if needed
