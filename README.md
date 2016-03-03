shkoma
======
Protein data analysis.


## Requirements
* Python 3 with installed `numpy`, `biopython` and `rpy2` packages
* R with installed `Peptides` package

## Usage example:

```python
from shkoma import correlation

# 1. load main experiment data
main_data = correlation.load_main_data_from_csv('experiment-28.csv')

# 2. construct list of proteins using main_data
proteins = correlation.construct_proteins(main_data)
correlation.fill_protein_sequences(proteins)
correlation.save_proteins_to_csv(proteins, 'proteins.csv')

# 3. alternate possibility: load proteins from .csv
proteins = correlation.load_proteins_from_csv('proteins.csv')

# 4. construct list of protein records
protein_records = correlation.construct_protein_records(proteins, main_data)

# 5. fill missed peptide records and computational protein parameters (for each protein)
correlation.fill_missed_peptide_records(protein_records)
correlation.fill_protein_parameters(protein_records)

# 6. save results to folder (separate file for each protein record)
correlation.save_protein_records_to_folder('experiment-28')
```