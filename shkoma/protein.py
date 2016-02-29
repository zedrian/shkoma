class Protein:
    def __init__(self, id, name, sequence=None):
        self.id = id
        self.name = name
        self.sequence = sequence

    def __str__(self):
        return '  ID: ' + self.id + '\n' + \
               '  Name: ' + self.name + '\n' + \
               '  Sequence: ' + self.sequence + '\n'

    def __repr__(self):
        return 'Protein general information:\n' + Protein.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.id == other.id


def find_protein_with_id(proteins, id):
    for protein in proteins:
        if protein.id == id:
            return protein

    return None
