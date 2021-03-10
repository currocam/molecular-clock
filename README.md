# Molecular clock

This python script allows you to compare different DNA sequences and elaborate a distance matrix without using the specific libraries developed for this. It was created with educational purposes. It may have applications in phylogenetic studies.

## General info

The script work as follows:

- Load the DNA sequences from the .fasta file
- Elaborate a DataFrame that contains appropriate info
- Elaborate a distance matrix using a Levenshtein distance and a Hamming distance calculator (in order to optimize execution time)

## Technologies
- Python 3.8.5
- Jupyter Notebook
## Example of use
An example of use is available in the "CytB_Birds.ipynb" file. 