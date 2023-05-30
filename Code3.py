from Bio import SeqIO

# Read the query sequence from the FASTA file
query_sequence = SeqIO.read("gene_sequence.fasta", "fasta")

# Read the reference sequences (chromosomes) from the FASTA file
reference_sequences = SeqIO.to_dict(SeqIO.parse("reference_genome.fasta", "fasta"))

# Calculate pairwise identity for each reference sequence
pairwise_identities = []
for chromosome_id, chromosome_sequence in reference_sequences.items():
    length = len(query_sequence)
    matches = sum(a == b for a, b in zip(query_sequence, chromosome_sequence))
    identity = (matches / length) * 100
    pairwise_identities.append((chromosome_id, identity))

# Find the best matching chromosome based on pairwise identity
best_matching_chromosome = max(pairwise_identities, key=lambda x: x[1])

# Print the results
print("Pairwise Identity for each chromosome:")
for chromosome_id, identity in pairwise_identities:
    print(f"Chromosome {chromosome_id}: {identity:.2f}%")

print("Best matching chromosome:", best_matching_chromosome[0])
print("Pairwise Identity:", best_matching_chromosome[1], "%")
