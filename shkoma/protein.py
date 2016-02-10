class Protein:
    def __init__(self, id='wrong id', name='wrong name', mw=0.0, pI=0.0, M=0.0, Z=0.0, sequence=''):
        self.id = id
        self.name = name
        self.mw = mw
        self.pI = pI
        self.M = M
        self.Z = Z
        self.sequence = sequence

    def __str__(self):
        return '(id: \'' + self.id + '\', name: \'' + self.name + '\'' + ', mw: ' + str(self.mw) + ', pI: ' + str(self.pI) + ', M: ' + str(self.M) + ', Z: ' \
               + str(self.Z) + ', sequence: \'' + self.sequence + '\')'

    def __repr__(self):
        return 'protein ' + Protein.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.id == other.id
