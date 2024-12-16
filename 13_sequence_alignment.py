# Module 13: Sequence Alignments

"""
Learning Objectives:
- Understand the concept of pairwise sequence alignments.
- Learn to perform global and local alignments using Bio.Align.
- Explore scoring schemes and advanced alignment parameters.

Key Concepts Covered:
1. Pairwise Sequence Alignments
2. Global Alignment
3. Local Alignment
4. Custom Scoring Schemes
5. Application of Alignments in Bioinformatics
"""

# 1. Pairwise Sequence Alignments
# Demonstration of aligning two sequences using Bio.Align
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Sample DNA sequences
dna_seq1 = "ACTGTTAC"
dna_seq2 = "ACTGTACT"

# Perform global alignment
alignments = pairwise2.align.globalxx(dna_seq1, dna_seq2)
print("Global Alignments:")
for alignment in alignments:
    print(format_alignment(*alignment))

# Perform local alignment
alignments = pairwise2.align.localxx(dna_seq1, dna_seq2)
print("\nLocal Alignments:")
for alignment in alignments:
    print(format_alignment(*alignment))

# 2. Global and Local Alignment Techniques
# Bio.Align provides functions for both types of alignments:
# - Global: Align sequences end-to-end.
# - Local: Find regions of highest similarity.

# 3. Custom Scoring Schemes
# Adjusting scoring schemes for match, mismatch, and gap penalties.

# Scoring example with custom match/mismatch/gap scores
dna_seq3 = "ACGT"
dna_seq4 = "AGT"
custom_alignments = pairwise2.align.globalms(
    dna_seq3, dna_seq4, match=2, mismatch=-1, open=-2, extend=-0.5
)
print("\nCustom Scoring - Global Alignment:")
for alignment in custom_alignments:
    print(format_alignment(*alignment))

# Practical Examples with Protein Sequences
protein_seq1 = "MVLSPADKTNVKAAW"
protein_seq2 = "MVLSAADKTNVKAAW"
protein_alignments = pairwise2.align.globalxx(protein_seq1, protein_seq2)
print("\nProtein Sequence Alignment:")
for alignment in protein_alignments:
    print(format_alignment(*alignment))

# Practical Exercise Questions
"""
1. Pairwise Alignment Challenge:
   Q1: Align the sequences "AGTACGCA" and "TATGC" using global alignment.
       Print the aligned sequences with their scores.

2. Local Alignment Challenge:
   Q2: Perform a local alignment on "GATTACA" and "GCATGCU".
       Identify the region of highest similarity.

3. Custom Scoring Challenge:
   Q3: Use a custom scoring scheme with match = 3, mismatch = -2, 
       gap open penalty = -3, and gap extend penalty = -1 to align 
       "ACCGT" and "ACG". Print the alignment and score.

4. Protein Alignment Challenge:
   Q4: Align the protein sequences "MKVQA" and "MKLQA" using global 
       alignment. Interpret the alignment results.
"""


def pairwise_alignment(seq1, seq2):
    """
    Performs global and local pairwise sequence alignments.

    Args:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence

    Returns:
        None: Prints alignments
    """
    # Perform global alignment
    print("Global Alignment:")
    global_alignments = pairwise2.align.globalxx(seq1, seq2)
    for alignment in global_alignments:
        print(format_alignment(*alignment))

    # Perform local alignment
    print("\nLocal Alignment:")
    local_alignments = pairwise2.align.localxx(seq1, seq2)
    for alignment in local_alignments:
        print(format_alignment(*alignment))


def custom_scoring_alignment(seq1, seq2, match, mismatch, open_gap, extend_gap):
    """
    Performs alignment using a custom scoring scheme.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        match (int): Score for matches
        mismatch (int): Penalty for mismatches
        open_gap (int): Gap opening penalty
        extend_gap (float): Gap extension penalty

    Returns:
        None: Prints alignments
    """
    print("\nCustom Scoring Alignment:")
    alignments = pairwise2.align.globalms(
        seq1, seq2, match, mismatch, open_gap, extend_gap
    )
    for alignment in alignments:
        print(format_alignment(*alignment))


def alignment_challenges():
    """
    Runs the exercise challenges for sequence alignment.
    """
    # Challenge 1: Global Alignment
    print("Challenge 1: Global Alignment")
    pairwise_alignment("AGTACGCA", "TATGC")

    # Challenge 2: Local Alignment
    print("\nChallenge 2: Local Alignment")
    pairwise_alignment("GATTACA", "GCATGCU")

    # Challenge 3: Custom Scoring
    print("\nChallenge 3: Custom Scoring")
    custom_scoring_alignment(
        "ACCGT", "ACG", match=3, mismatch=-2, open_gap=-3, extend_gap=-1
    )

    # Challenge 4: Protein Sequence Alignment
    print("\nChallenge 4: Protein Sequence Alignment")
    pairwise_alignment("MKVQA", "MKLQA")


# Uncomment the following lines to test the module
# alignment_challenges()
