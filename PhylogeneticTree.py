import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from Bio import SeqIO
import pandas as pd
from reportlab.lib import colors
from reportlab.platypus import Table, TableStyle
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas


def validate_dna(seq):
    """
    Validates if a sequence contains only valid DNA nucleotides (A, C, G, T).

    Arguments:
        seq (str): DNA sequence string.

    Returns:
        bool: True if valid DNA sequence, False otherwise.
    """
    return all(c in 'ACGTacgt' for c in seq)


def load_sequences(file=None, seqs=None):
    """
    Loads DNA sequences either from a FASTA file or a list of sequences,
    validating each sequence contains only valid DNA nucleotides.

    Arguments:
        file (str, optional): Path to FASTA file containing sequences.
        seqs (list of str, optional): List of sequence strings.

    Returns:
        tuple: (sequences, labels)
            sequences (list of str): Validated DNA sequences in uppercase.
            labels (list of str): Corresponding sequence labels or IDs.

    Raises:
        ValueError: If fewer than two valid sequences are provided,
                    or if neither file nor sequences are given.
    """

    sequences = []
    labels = []
    if file:
        for record in SeqIO.parse(file, "fasta"):
            seq = str(record.seq).upper()
            if validate_dna(seq):
                sequences.append(seq)
                labels.append(record.id)
        if len(sequences) < 2:
            raise ValueError("At least two valid DNA sequences required.")
    elif seqs:
        sequences = [s.upper() for s in seqs if validate_dna(s)]
        labels = [f"Seq{i+1}" for i in range(len(sequences))]
        if len(sequences) < 2:
            raise ValueError("At least two valid DNA sequences required.")
    else:
        raise ValueError("Provide either --file or --seqs input.")
    return sequences, labels


def compute_distance(seq1, seq2, match=1, mismatch=0, gap=-1):
    """
    Computes a normalized distance between two DNA sequences based on
    a global alignment scoring with match, mismatch, and gap penalties.

    Arguments:
        seq1 (str): First DNA sequence.
        seq2 (str): Second DNA sequence.
        match (int, optional): Score for a match (default: 1).
        mismatch (int, optional): Penalty for mismatch (default: 0).
        gap (int, optional): Penalty for gap (default: -1).

    Returns:
        float: Normalized distance between sequences (0 to 1),
               where 0 means identical and 1 means completely different.
    """

    len1, len2 = len(seq1), len(seq2)
    dp = np.zeros((len1 + 1, len2 + 1))

    for i in range(len1 + 1):
        dp[i][0] = i * gap
    for j in range(len2 + 1):
        dp[0][j] = j * gap

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            diag = dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = dp[i - 1][j] + gap
            insert = dp[i][j - 1] + gap
            dp[i][j] = max(diag, delete, insert)

    max_score = min(len1, len2) * match
    identity = dp[len1][len2] / max_score if max_score != 0 else 0
    return 1 - identity


def build_distance_matrix(sequences, match=1, mismatch=0, gap=-1):
    """
    Builds a symmetric pairwise distance matrix for a list of sequences
    using the compute_distance function.

    Arguments:
        sequences (list of str): List of DNA sequences.
        match (int, optional): Match score (default: 1).
        mismatch (int, optional): Mismatch penalty (default: 0).
        gap (int, optional): Gap penalty (default: -1).

    Returns:
        numpy.ndarray: 2D symmetric matrix of distances.
    """

    n = len(sequences)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist = compute_distance(sequences[i], sequences[j], match, mismatch, gap)
            matrix[i][j] = matrix[j][i] = dist
    return matrix


def upgma_clustering(dist_matrix, labels, output_img="upgma_tree.png"):
    """
    Performs UPGMA hierarchical clustering on a distance matrix,
    generates and saves a dendrogram plot and writes clustering results.

    Arguments:
        dist_matrix (numpy.ndarray): Symmetric pairwise distance matrix.
        labels (list of str): Sequence labels corresponding to matrix rows.
        output_img (str, optional): Filename for saving dendrogram image.
    """

    condensed = squareform(dist_matrix)
    linkage_matrix = linkage(condensed, method='average')

    plt.figure(figsize=(12, 6))
    dendrogram(linkage_matrix, labels=labels, orientation='top')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.title("UPGMA Phylogenetic Tree")

    # Removes y-axis ticks and label for a cleaner look
    plt.gca().yaxis.set_ticks([])
    plt.gca().set_ylabel('')

    plt.tight_layout()
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    plt.close()

    with open("upgma_results.txt", "w") as f:
        f.write("Input Sequences/Distances:\n")
        for i, label in enumerate(labels):
            f.write(f"{label}: {dist_matrix[i]}\n")
        f.write("\nUPGMA Clustering Results:\n")
        np.savetxt(f, linkage_matrix, fmt='%1.4f')



def generate_pdf_report(sequences, labels, dist_matrix, tree_img="upgma_tree.png", pdf_file="phylo_report.pdf"):
    """
    Generates a PDF report including input sequences, pairwise distance matrix,
    and UPGMA dendrogram tree image, ensuring proper layout across pages.

    Arguments:
        sequences (list of str): List of input DNA sequences (optional).
        labels (list of str): Corresponding sequence labels.
        dist_matrix (numpy.ndarray): Pairwise distance matrix.
        tree_img (str, optional): Path to UPGMA dendrogram image file.
        pdf_file (str, optional): Output PDF filename.
    """

    c = canvas.Canvas(pdf_file, pagesize=letter)
    width, height = letter
    margin = 50
    y = height - margin

    # Title
    c.setFont("Helvetica-Bold", 16)
    c.drawString(margin, y, "UPGMA Phylogenetic Analysis Report")
    y -= 40

    if sequences:
        # Input sequences section
        c.setFont("Helvetica-Bold", 12)
        c.drawString(margin, y, "Input Sequences:")
        y -= 30
        c.setFont("Helvetica", 9)

        for i, seq in enumerate(sequences):
            line = f"{labels[i]}: {seq}"
            max_chars_per_line = 90
            while len(line) > 0:
                c.drawString(margin, y, line[:max_chars_per_line])
                line = line[max_chars_per_line:]
                y -= 12
                if y < 100:
                    c.showPage()
                    y = height - margin
                    c.setFont("Helvetica", 9)

            y -= 6  # spacing between sequences

        c.showPage()
        y = height - margin
    else:
        # Matrix was preloaded; no sequence input
        c.setFont("Helvetica", 11)
        c.drawString(margin, y, "Note: Sequences were not provided directly; using a distance matrix given a priori.")
        y -= 30

    # Distance matrix section
    c.setFont("Helvetica-Bold", 14)
    c.drawString(margin, y, "Pairwise Distance Matrix:")
    y -= 30

    short_labels = [lbl if len(lbl) <= 12 else lbl[:9] + "..." for lbl in labels]
    data = [[""] + short_labels]
    for i, row_label in enumerate(short_labels):
        row_vals = [f"{v:.2f}" for v in dist_matrix[i]]
        data.append([row_label] + row_vals)

    table = Table(data, repeatRows=1)
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 8),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 6),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.whitesmoke, colors.lightgrey])
    ])
    table.setStyle(style)

    max_table_width = width - 2 * margin
    table_width, table_height = table.wrap(0, 0)
    scale_factor = min(1, max_table_width / table_width)
    num_cols = len(data[0])
    col_widths = [(table_width / num_cols) * scale_factor for _ in range(num_cols)]

    table = Table(data, colWidths=col_widths)
    table.setStyle(style)
    _, table_height = table.wrapOn(c, max_table_width, height)
    table.drawOn(c, margin, y - table_height)

    # Dendrogram page
    c.showPage()
    c.setFont("Helvetica-Bold", 14)
    c.drawString(margin, height - margin, "UPGMA Tree:")
    try:
        c.drawImage(tree_img, margin, height - margin - 320, width=500, height=300)
    except Exception as e:
        c.drawString(margin, height - margin - 20, f"Error loading tree image: {e}")

    c.save()
    print(f"PDF report saved to: {pdf_file}")


def main():
    """
    Parses command-line arguments to load sequences or distance matrix,
    computes distance matrix if needed, performs clustering, and generates reports.
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
        df = pd.read_csv(args.matrix, index_col=0)
        labels = df.index.tolist()
        dist_matrix = df.values
        if dist_matrix.shape[0] != dist_matrix.shape[1]:
            raise ValueError("Distance matrix must be square.")
        if not np.allclose(dist_matrix, dist_matrix.T):
            raise ValueError("Distance matrix must be symmetric.")
        print("Loaded distance matrix from CSV.")
        sequences = None
    else:
        sequences, labels = load_sequences(file=args.file, seqs=args.seqs)
        dist_matrix = build_distance_matrix(sequences, args.match, args.mismatch, args.gap)
        pd.DataFrame(dist_matrix, index=labels, columns=labels).to_csv("distance_matrix.csv")
        print("Computed pairwise distances and saved to distance_matrix.csv.")

    upgma_clustering(dist_matrix, labels)
    print("Generated dendrogram and clustering results.")

    generate_pdf_report(
        sequences=sequences if sequences else None,
        labels=labels,
        dist_matrix=dist_matrix
    )


if __name__ == "__main__":
    main()
