# UPGMA Phylogenetic Tree Builder

This project allows you to construct a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm. It can compute distances from a set of DNA sequences or use a precomputed distance matrix. The resulting tree is visualized, and the analysis is saved in a report.

## Features

- Flexible Input: Load DNA sequences from a FASTA file, input them manually, or provide a custom distance matrix via CSV.
- Validation of DNA sequences: Ensures that sequences contain only valid DNA characters (A, C, G, T).
- Distance Calculation: Computes pairwise genetic distances using Needleman-Wunsch global alignment with customizable scoring parameters.
- UPGMA Clustering: Builds a phylogenetic tree based on genetic distances.
- Dendrogram Visualization: Generates and saves a dendrogram tree image representing phylogenetic relationships.
- PDF Report Generation: Outputs a styled PDF containing the input sequences, pairwise distance matrix and the dendrogram image.
- Result Export: Saves the computed distance matrix, dendrogram image, and clustering linkage data for further use.
  
## How to Run

### Step 1: Prepare Your Environment
Make sure you have the required Python packages installed:
```bash
pip install biopython matplotlib scipy pandas
```
### Step 2: Run the Program
Run the script from the command line:
```bash
python PhylogeneticTree.py --file sequences.fasta
```
Or with direct sequence input:

```bash
python PhylogeneticTree.py --seqs ACTGCT AGTGCT AGTGAT
```
Or by using a custom distance matrix:
```bash
python PhylogeneticTree.py --matrix distances.csv
```
### Step 3: View Outputs
After running, the following files will be generated:

- upgma_tree.png: A dendrogram representing the phylogenetic tree.
- upgma_results.txt: Text file with linkage data and input matrix/labels.
- distance_matrix.csv (only for sequence input): The computed pairwise distances.
- phylo_report.pdf: A comprehensive PDF report with sequences, distance matrix, and tree.
  
## Input Formats

### FASTA File
A standard FASTA file with DNA sequences:
```bash
>Seq1
ACTGCTAG
>Seq2
ACTGCTAC
>Seq3
ACTGATAC
```

### Distance Matrix CSV
A square symmetric CSV file:
```bash
,Seq1,Seq2,Seq3
Seq1,0,1,2
Seq2,1,0,3
Seq3,2,3,0
```

## Customization

You can adjust scoring for alignment-based distance calculation using:
```bash
--match 1 --mismatch 0 --gap -1
```
These options apply when input sequences are used (not when using a matrix).

## Example
```bash
python PhylogeneticTree.py --seqs ACTG ACTA ACGT --match 1 --mismatch 0 --gap -1
```
This will compute distances, build a UPGMA tree, and save all outputs.
