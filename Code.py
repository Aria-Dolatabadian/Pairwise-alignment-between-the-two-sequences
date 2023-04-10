from Bio import pairwise2
from Bio import SeqIO
import csv

# Load the two DNA sequences in fasta format
seq1 = SeqIO.read('sequence1.fasta', 'fasta')
seq2 = SeqIO.read('sequence2.fasta', 'fasta')

# Calculate the pairwise alignment between the two sequences
alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)

# Extract the alignment information and write to a CSV file
with open('alignments.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Sequence 1 ID', 'Sequence 1 Alignment', 'Sequence 2 ID', 'Sequence 2 Alignment', 'Alignment Score'])
    for alignment in alignments:
        writer.writerow([seq1.id, alignment.seqA, seq2.id, alignment.seqB, alignment.score])

# Get the best alignment score and alignment strings
best_alignment = max(alignments, key=lambda x: x.score)
score = best_alignment.score
alignment_seq1 = best_alignment.seqA
alignment_seq2 = best_alignment.seqB

# Calculate the positions of the aligned nucleotides
aligned_positions = [(i, j) for i, j in zip(range(len(alignment_seq1)), range(len(alignment_seq2))) if alignment_seq1[i] != '-' and alignment_seq2[j] != '-']

# Plot the positions of the aligned nucleotides
import matplotlib.pyplot as plt

x = [i for i, j in aligned_positions]
y = [j for i, j in aligned_positions]
plt.scatter(x, y)
plt.xlabel('Sequence 1')
plt.ylabel('Sequence 2')
plt.title(f'Pairwise Alignment Score: {score}')
plt.show()








