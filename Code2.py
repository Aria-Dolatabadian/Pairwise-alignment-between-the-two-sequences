from Bio import pairwise2
from Bio import SeqIO
import matplotlib.pyplot as plt

# Define colors for highlighting aligned nucleotides
match_color = '#00FF00'
mismatch_color = 'red'
gap_color = '#808080'


# Define a function to plot the two sequences with aligned nucleotides highlighted
def plot_alignment(seq1, seq2, aligned_positions):
    fig, ax = plt.subplots(figsize=(15, 5))

    # Plot sequence 1
    ax.text(0, 0.5, seq1, fontfamily='monospace', fontsize=14)

    # Plot sequence 2
    ax.text(0, -0.5, seq2, fontfamily='monospace', fontsize=14)

    # Highlight aligned nucleotides
    for i, j in aligned_positions:
        if seq1[i] == seq2[j]:
            color = match_color
        elif seq1[i] == '-' or seq2[j] == '-':
            color = gap_color
        else:
            color = mismatch_color
        ax.add_patch(plt.Rectangle((i - 0.5, -0.6), 1, 1.2, color=color, alpha=0.5))

    # Set axis limits and remove tick marks
    ax.set_xlim([-0.5, len(seq1) - 0.5])
    ax.set_ylim([-1, 1])
    ax.set_xticks([])
    ax.set_yticks([])

    plt.show()


# Load the two DNA sequences in fasta format
seq1 = SeqIO.read('sequence1.fasta', 'fasta')
seq2 = SeqIO.read('sequence2.fasta', 'fasta')

# Calculate the pairwise alignment between the two sequences
alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)

# Get the best alignment score and alignment strings
best_alignment = max(alignments, key=lambda x: x.score)
score = best_alignment.score
alignment_seq1 = best_alignment.seqA
alignment_seq2 = best_alignment.seqB

# Calculate the positions of the aligned nucleotides
aligned_positions = [(i, j) for i, j in zip(range(len(alignment_seq1)), range(len(alignment_seq2))) if
                     alignment_seq1[i] != '-' and alignment_seq2[j] != '-']

# Call the function to plot the alignment
plot_alignment(alignment_seq1, alignment_seq2, aligned_positions)
