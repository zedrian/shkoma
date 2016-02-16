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
        return '(analysis_name: \'' + self.analysis_name + '\', score: ' + str(self.score) + \
               ', reverse_score: ' + str(self.reverse_score) + \
               ', percent_of_scored_peak_intensity: ' + str(self.percent_of_scored_peak_intensity) + \
               ', total_intensity: ' + str(self.total_intensity) + \
               ', precursor_averagine_chi_squared: ' + str(self.precursor_averagine_chi_squared) + \
               ', retention_time_min: ' + str(self.retention_time_min) + \
               ', chromatographic_peak_width_in_seconds: ' + str(self.chromatographic_peak_width_in_seconds) + ')'

    def __repr__(self):
        return PeptideMatch.__str__(self)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.analysis_name == other.analysis_name  # TODO: add other checks if needed
