import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from Bio import SeqIO
import pandas as pd


def validate_dna(seq):
    """
    Validates whether a DNA sequence contains only valid nucleotides (A, C, G, T).

    Arguments:
        seq (str): DNA sequence to validate.

    Returns:
        bool: True if sequence is valid, False otherwise.
    """
    return all(c in 'ACGTacgt' for c in seq)


def load_sequences(file=None, seqs=None):
    """
    Loads DNA sequences from a FASTA file or a list of input strings.

    Arguments:
        file (str, optional): Path to a FASTA file containing DNA sequences.
        seqs (list[str], optional): List of DNA sequence strings.

    Returns:
        list[str]: List of validated DNA sequences in uppercase.

    Raises:
        ValueError: If fewer than two valid sequences are found or no input is provided.
    """
    sequences = []
    if file:
        for record in SeqIO.parse(file, "fasta"):
            seq = str(record.seq).upper()
            if validate_dna(seq):
                sequences.append(seq)
        if len(sequences) < 2:
            raise ValueError("At least two valid DNA sequences required.")
    elif seqs:
        sequences = [s.upper() for s in seqs if validate_dna(s)]
        if len(sequences) < 2:
            raise ValueError("At least two valid DNA sequences required.")
    else:
        raise ValueError("Provide either --file or --seqs input.")
    return sequences


def compute_distance(seq1, seq2, match=1, mismatch=0, gap=-1):
    """
    Computes a simple distance between two sequences based on character matches.

    Arguments:
        seq1 (str): First DNA sequence.
        seq2 (str): Second DNA sequence.
        match (int): Score for matching characters.
        mismatch (int): Penalty for mismatched characters.
        gap (int): Penalty for gaps (currently unused in this function).

    Returns:
        float: Normalized distance score (0 = identical, 1 = completely different).
    """
    matches = sum(a == b for a, b in zip(seq1, seq2))
    length = max(len(seq1), len(seq2))
    return 1 - (matches / length)


def build_distance_matrix(sequences, match=1, mismatch=0, gap=-1):
    """
    Constructs a symmetric pairwise distance matrix for a list of sequences.

    Arguments:
        sequences (list[str]): List of DNA sequences.
        match (int): Score for matching characters.
        mismatch (int): Penalty for mismatched characters.
        gap (int): Penalty for gaps (currently unused in this function).

    Returns:
        ndarray: Symmetric NxN distance matrix where N is the number of sequences.
    """
    n = len(sequences)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist = compute_distance(sequences[i], sequences[j], match, mismatch, gap)
            matrix[i][j] = matrix[j][i] = dist
    return matrix


def upgma_clustering(dist_matrix, labels):
    """
    Performs UPGMA hierarchical clustering on a distance matrix and saves a dendrogram plot.

    Arguments:
        dist_matrix (ndarray): NxN symmetric distance matrix.
        labels (list[str]): List of labels corresponding to the sequences.

    Saves:
        upgma_tree.png: Dendrogram plot of the clustering.
        upgma_results.txt: Linkage matrix and input matrix/labels.
    """
    condensed = squareform(dist_matrix)
    linkage_matrix = linkage(condensed, method='average')

    plt.figure(figsize=(10, 5))
    dendrogram(linkage_matrix, labels=labels, orientation='top')
    plt.title("UPGMA Phylogenetic Tree")
    plt.ylabel("Genetic Distance")
    plt.savefig("upgma_tree.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save results to file
    with open("upgma_results.txt", "w") as f:
        f.write("Input Sequences/Distances:\n")
        if len(labels) == len(dist_matrix):
            for i, label in enumerate(labels):
                f.write(f"{label}: {dist_matrix[i]}\n")
        f.write("\nUPGMA Clustering Results:\n")
        np.savetxt(f, linkage_matrix, fmt='%1.4f')


def main():
    """
    Main function to parse arguments, process input data, perform UPGMA clustering, and save results.

    Input options (mutually exclusive):
        --file (str): Path to FASTA file containing sequences.
        --seqs (list[str]): Space-separated DNA sequences.
        --matrix (str): Path to CSV file containing a distance matrix.

    Optional arguments:
        --match (int): Match score (default: 1).
        --mismatch (int): Mismatch penalty (default: 0).
        --gap (int): Gap penalty (default: -1).

    Outputs:
        - upgma_tree.png: Image of the dendrogram.
        - upgma_results.txt: Clustering results.
        - distance_matrix.csv: Computed pairwise distances (if sequences are used).
    """
    parser = argparse.ArgumentParser(description="UPGMA Phylogenetic Tree Builder")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--file", help="FASTA file containing sequences")
    input_group.add_argument("--seqs", nargs='+', help="Space-separated sequences")
    input_group.add_argument("--matrix", help="CSV file with distance matrix")

    parser.add_argument("--match", type=int, default=1, help="Match score (default: 1)")
    parser.add_argument("--mismatch", type=int, default=0, help="Mismatch penalty (default: 0)")
    parser.add_argument("--gap", type=int, default=-1, help="Gap penalty (default: -1)")

    args = parser.parse_args()

    if args.matrix:
        # Load matrix and ensure it's square
        df = pd.read_csv(args.matrix, index_col=0)
        labels = df.index.tolist()
        dist_matrix = df.values

        if dist_matrix.shape[0] != dist_matrix.shape[1]:
            raise ValueError("Distance matrix must be square.")

        if not np.allclose(dist_matrix, dist_matrix.T):
            raise ValueError("Distance matrix must be symmetric.")

        print("Loaded distance matrix from CSV.")

    else:
        # Load and validate sequences
        sequences = load_sequences(file=args.file, seqs=args.seqs)
        labels = [f"Seq{i + 1}" for i in range(len(sequences))]
        dist_matrix = build_distance_matrix(sequences, args.match, args.mismatch, args.gap)
        pd.DataFrame(dist_matrix, index=labels, columns=labels).to_csv("distance_matrix.csv")
        print("Computed pairwise distances and saved to distance_matrix.csv.")

    upgma_clustering(dist_matrix, labels)
    print("Analysis complete. Results saved to:")
    print("- upgma_results.txt (clustering data)")
    print("- upgma_tree.png (dendrogram)")
    if not args.matrix:
        print("- distance_matrix.csv (pairwise distances)")



if __name__ == "__main__":
    main()