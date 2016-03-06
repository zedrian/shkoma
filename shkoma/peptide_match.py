class PeptideMatch:
    def __init__(self, analysis_name='', score=0.0, reverse_score=0.0,
                 percent_of_scored_peak_intensity=0.0, total_intensity=0.0,
                 precursor_averagine_chi_squared=0.0, retention_time_min=0.0,
                 chromatographic_peak_width_in_seconds=0.0):
        self.analysis_name = analysis_name
        self.score = score
        self.reverse_score = reverse_score
        self.percent_of_scored_peak_intensity = percent_of_scored_peak_intensity
        self.total_intensity = total_intensity
        self.precursor_averagine_chi_squared = precursor_averagine_chi_squared
        self.retention_time_min = retention_time_min
        self.chromatographic_peak_width_in_seconds = chromatographic_peak_width_in_seconds

    def __str__(self):
        return '    Analysis name: {0}\n'.format(self.analysis_name) + \
               '    Score: {0}\n'.format(self.score) + \
               '    Reverse score: {0}\n'.format(self.reverse_score) + \
               '    Percent of scored peak intensity: {0}\n'.format(self.percent_of_scored_peak_intensity) + \
               '    Total intensity: {0}\n'.format(self.total_intensity) + \
               '    Precursor averagine chi squared: {0}\n'.format(self.precursor_averagine_chi_squared) + \
               '    Retention time min: {0}\n'.format(self.retention_time_min) + \
               '    Chromatographic peak width (in seconds): {0}\n'.format(self.chromatographic_peak_width_in_seconds)

    def __repr__(self):
        return '\nPeptide match:\n' + PeptideMatch.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.analysis_name == other.analysis_name  # TODO: add other checks if needed
