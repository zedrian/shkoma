from shkoma.peptide import Peptide


class PeptideRecord:
    def __init__(self, peptide=Peptide, matches=[]):
        self.peptide = peptide
        self.matches = matches

    def __str__(self):
        return 'peptide: ' + str(self.peptide) + ', matches: ' + str(self.matches)

    def __repr__(self):
        return 'peptide record (' + PeptideRecord.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.peptide == other.peptide and \
               self.matches == other.matches
