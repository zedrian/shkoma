class ProteinParameters:
    def __init__(self, sequence):
        self.sequence = sequence
        self.sequence_length = len(sequence)
        self.M = 0
        self.Z = 0

    def __str__(self):
        return '  Sequence length: ' + str(self.sequence_length) + '\n' + \
               '  M: {0:.3f}\n'.format(self.M) + \
               '  Z: {0:.3f}\n'.format(self.Z)

    def __repr__(self):
        return 'Protein computational parameters:\n' + ProteinParameters.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.sequence == other.sequence
