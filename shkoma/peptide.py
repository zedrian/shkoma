class Peptide:
    def __init__(self, sequence='', pI=0.0, aminoacid_compound=None,
                 logP=0.0, aminoacid_matrix=None, aminoacid_types=None,
                 sequence_place=None, structure_place=None):
        self.sequence = sequence
        self.pI = pI
        self.aminoacid_compound = aminoacid_compound
        self.logP = logP
        self.aminoacid_matrix = aminoacid_matrix
        self.aminoacid_types = aminoacid_types
        self.sequence_place = sequence_place
        self.structure_place = structure_place

    def __str__(self):
        return '(sequence: \'' + self.sequence + '\', pI: ' + str(self.pI) + \
               ', aminoacid_compound: ' + str(self.aminoacid_compound) + \
               ', logP: ' + str(self.logP) + ', aminoacid_matrix: ' + str(self.aminoacid_matrix) + \
               ', aminoacid_types: ' + str(self.aminoacid_types) + \
               ', sequence_place: ' + str(self.sequence_place) + \
               ', structure_place: ' + str(self.structure_place) + ')'

    def __repr__(self):
        return 'peptide ' + Peptide.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.sequence == other.sequence  # TODO: add other checks if needed
