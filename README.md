# UPGMA Phylogenetic Tree Builder

This project allows you to construct a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm. It can compute distances from a set of DNA sequences or use a precomputed distance matrix. The resulting tree is visualized and saved as an image, and results are logged in a text file.

## Features

- Flexible Input: Load DNA sequences from a FASTA file, input them manually, or provide a custom distance matrix via CSV.
- Validation of DNA sequences: Ensures that sequences contain only valid DNA characters (A, C, G, T).
- Distance Calculation: Automatically computes a pairwise distance matrix from input sequences.
- UPGMA Clustering: Builds a phylogenetic tree based on genetic distances.
- Dendrogram Visualization: Saves the tree structure as a high-resolution PNG image.
- Result Export: Saves a detailed clustering report and (if applicable) the computed distance matrix.
  
## How to Run

### Step 1: Prepare Your Environment
Make sure you have the required Python packages installed:
```bash
pip install biopython matplotlib scipy pandas
```
### Step 2: Run the Program
Run the script from the command line:
```bash
python PhyloTree.py --file sequences.fasta
```
Or with direct sequence input:

```bash
python PhyloTree.py --seqs ACTGCT AGTGCT AGTGAT
```
Or by using a custom distance matrix:
```bash
python PhyloTree.py --matrix distances.csv
```
### Step 3: View Outputs
After running, the following files will be generated:

- upgma_tree.png: A dendrogram representing the phylogenetic tree.
- upgma_results.txt: Text file with linkage data and input matrix/labels.
- distance_matrix.csv (only for sequence input): The computed pairwise distances.
  
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
Seq1,0,0.1,0.2
Seq2,0.1,0,0.3
Seq3,0.2,0.3,0
```

## Customization

You can adjust scoring for alignment-based distance calculation using:
```bash
--match 1 --mismatch 0 --gap -1
```
These options apply when input sequences are used (not when using a matrix).

## Example
```bash
python PhyloTree.py --seqs ACTG ACTA ACGT --match 1 --mismatch 0 --gap -1
```
This will compute distances, build a UPGMA tree, and save all outputs to disk.
