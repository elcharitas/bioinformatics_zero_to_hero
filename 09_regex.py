# Module 9: Regular Expressions

"""
Learning Objectives:
- Understand the basics of the `re` module for pattern matching.
- Learn how to use regular expressions to extract motifs from DNA sequences.

Key Concepts Covered:
1. Regular Expression Basics
2. Using the `re` Module
3. Extracting Biological Motifs
"""

import re

# 1. Regular Expression Basics
# Regular expressions (regex) are used to match patterns in text.
# Example: A simple regex to match the pattern 'ATG' in a DNA sequence.
dna_sequence = "ATGCGTACGTTGATGCTAG"
pattern = "ATG"

# Using `re.findall` to find all occurrences of the pattern in the sequence
matches = re.findall(pattern, dna_sequence)
print(f"Matches for '{pattern}' in the DNA sequence: {matches}")

# Additional Example: Matching multiple patterns
# Finds both 'ATG' and 'GCT' in the sequence
pattern_multi = r"ATG|GCT"  # Using '|' to match either pattern
matches_multi = re.findall(pattern_multi, dna_sequence)
print(f"Matches for '{pattern_multi}' in the DNA sequence: {matches_multi}")

# 2. Using the `re` Module for Complex Patterns
# Regex to find sequences starting with 'ATG' and followed by any 3 nucleotides
dna_sequence = "ATGCGTACGTATGTTAGCTGATGACG"
pattern = r"ATG..."  # '...' matches any three characters
matches = re.findall(pattern, dna_sequence)
print(f"Matches for pattern '{pattern}': {matches}")

# Additional Example: Finding overlapping patterns
# Matches overlapping 'ATG's in the sequence
overlap_pattern = r"(?=(ATG))"  # Lookahead to allow overlaps
matches_overlap = re.findall(overlap_pattern, dna_sequence)
print(f"Overlapping matches for '{overlap_pattern}': {matches_overlap}")

# 3. Extracting Biological Motifs
# Example: Extracting a specific DNA motif like a promoter or a binding site.
dna_sequence = "TATAATCGTAGCATGTACGTGCGTAC"
pattern = r"TATA[ATGC]{4}"  # Matches 'TATA' followed by any 4 nucleotides
matches = re.findall(pattern, dna_sequence)
print(f"Biological motifs matching '{pattern}': {matches}")

# Additional Example: Extracting start codons and surrounding sequences
# Matches 'ATG' followed by exactly 5 nucleotides
dna_sequence = "ATGCGTACGTATGTTAGCTGATGACG"
pattern_start_codon = r"ATG[ATGC]{5}"
matches_start_codon = re.findall(pattern_start_codon, dna_sequence)
print(f"Start codons with surrounding sequence: {matches_start_codon}")

# Practical Exercise Questions
"""
1. Pattern Matching Challenge:
   Q1: Write a regex pattern to find all occurrences of 'CG' in a DNA sequence.
   Example Input: "CGTAGCGCGT"
   Expected Output: ['CG', 'CG', 'CG']

2. Motif Extraction Challenge:
   Q2: Write a regex pattern to identify sequences that start with 'A', end with 'T',
       and have exactly 5 characters in between.
   Example Input: "ATGCGTACGT"
   Expected Output: ['ATGCGT']

3. DNA Sequence Validation:
   Q3: Write a regex pattern to check if a sequence contains only valid DNA bases (A, T, C, G).

4. Overlapping Pattern Matching:
   Q4: Write a regex pattern to find all overlapping occurrences of 'ATG' in a DNA sequence.
"""


def find_pattern(dna_sequence, pattern):
    """
    Finds all matches of a regex pattern in a DNA sequence.

    Args:
        dna_sequence (str): The DNA sequence to search.
        pattern (str): The regex pattern to match.

    Returns:
        list: All matches of the pattern.
    """
    matches = re.findall(pattern, dna_sequence)
    return matches


def validate_dna_sequence(dna_sequence):
    """
    Validates if a DNA sequence contains only valid bases (A, T, C, G).

    Args:
        dna_sequence (str): The DNA sequence to validate.

    Returns:
        bool: True if valid, False otherwise.
    """
    pattern = r"^[ATCG]+$"  # Matches sequences containing only A, T, C, G
    return bool(re.match(pattern, dna_sequence))


def find_overlapping_patterns(dna_sequence, pattern):
    """
    Finds overlapping matches of a regex pattern in a DNA sequence.

    Args:
        dna_sequence (str): The DNA sequence to search.
        pattern (str): The regex pattern with lookahead for overlaps.

    Returns:
        list: All overlapping matches of the pattern.
    """
    matches = re.findall(pattern, dna_sequence)
    return matches


# Uncomment the following lines to test the functions
# print(find_pattern("CGTAGCGCGT", "CG"))
# print(validate_dna_sequence("ATGCGTACGTT"))
# print(find_overlapping_patterns("ATGATGATG", r"(?=(ATG))"))
