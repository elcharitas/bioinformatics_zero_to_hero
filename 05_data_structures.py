# Module 5: Lists, Tuples, and Dictionaries

"""
Learning Objectives:
- Understand the structure and use of lists, tuples, and dictionaries.
- Learn common operations and methods for these data structures.
- Explore applications in storing and manipulating biological sequence data.

Key Concepts Covered:
1. Lists
2. Tuples
3. Dictionaries
4. Applications in Bioinformatics
"""

# 1. Lists
# Lists are ordered, mutable collections of items.
# Example: Representing a DNA sequence as a list of bases.
dna_sequence = ["A", "T", "G", "C", "A"]  # List of DNA bases
print(f"Original DNA Sequence: {dna_sequence}")

# Modifying a list
# Adding a base to the end of the sequence
rna_sequence = dna_sequence.copy()  # Creating an RNA sequence for modification
rna_sequence[-1] = "U"  # Changing the last base for RNA representation
print(f"RNA Sequence: {rna_sequence}")

# List Methods
# Count the occurrences of a base in the sequence
print(f"Occurrences of 'A': {dna_sequence.count('A')}")
# Reverse the sequence
rna_sequence.reverse()
print(f"Reversed RNA Sequence: {rna_sequence}")

# Extending the list
additional_bases = ["G", "A"]
dna_sequence.extend(additional_bases)
print(f"Extended DNA Sequence: {dna_sequence}")

# Sorting the list
dna_sequence.sort()
print(f"Sorted DNA Sequence: {dna_sequence}")

# 2. Tuples
# Tuples are immutable, ordered collections.
# Example: Storing a codon (triplet of bases).
codon = ("A", "T", "G")
print(f"Codon: {codon}")

# Accessing elements in a tuple
print(f"First base in codon: {codon[0]}")

# Iterating over a tuple
for base in codon:
    print(f"Base in codon: {base}")

# Nested tuples
nested_codon = ("A", ("T", "G"))
print(f"Nested Codon: {nested_codon}")

# 3. Dictionaries
# Dictionaries store key-value pairs.
# Example: Mapping codons to their corresponding amino acids.
codon_to_amino_acid = {
    "ATG": "Methionine",
    "TGG": "Tryptophan",
    "TAA": "STOP",
    "TAG": "STOP",
    "TGA": "STOP",
}

print(f"Codon to Amino Acid Mapping: {codon_to_amino_acid}")

# Accessing a value using its key
print(f"Amino acid for ATG: {codon_to_amino_acid['ATG']}")

# Adding a new codon to the dictionary
codon_to_amino_acid["GGC"] = "Glycine"
print(f"Updated Codon Mapping: {codon_to_amino_acid}")

# Iterating over a dictionary
for codon, amino_acid in codon_to_amino_acid.items():
    print(f"Codon: {codon}, Amino Acid: {amino_acid}")

# Checking if a key exists
if "TAA" in codon_to_amino_acid:
    print("TAA is a valid stop codon.")

# Removing a key-value pair
codon_to_amino_acid.pop("TAG")
print(f"Updated Codon Mapping after removal: {codon_to_amino_acid}")

# Practical Applications
# Use case: Translating a DNA sequence into a protein sequence
protein_sequence = []
for i in range(0, len(dna_sequence) - 2, 3):  # Iterate through triplets
    codon = "".join(dna_sequence[i : i + 3])  # Join triplet bases into a codon
    protein_sequence.append(codon_to_amino_acid.get(codon, "Unknown"))

print(f"Protein Sequence: {protein_sequence}")

# Practical Exercise Questions
"""
1. List Manipulation:
   Q1: Write a script that takes a list of RNA bases and replaces all instances of "U" with "T" to convert it back to DNA.

2. Tuple and Dictionary Challenge:
   Q2: Create a tuple representing a codon and a dictionary mapping the codon to its amino acid.
       Print the tuple and retrieve its corresponding amino acid from the dictionary.

3. Advanced List Operations:
   Q3: Write a script to find the most frequent base in a DNA sequence.

4. Bioinformatics Application:
   Q4: Write a function that translates a full DNA sequence into a protein sequence using a codon-to-amino-acid dictionary.
"""


def rna_to_dna_conversion(rna_sequence):
    """
    Converts RNA sequence (list) to DNA sequence by replacing "U" with "T".

    Args:
        rna_sequence (list): List of RNA bases.
    Returns:
        list: DNA sequence.
    """
    return ["T" if base == "U" else base for base in rna_sequence]


def translate_dna_to_protein(dna_sequence, codon_map):
    """
    Translates a DNA sequence into a protein sequence using a codon-to-amino-acid dictionary.

    Args:
        dna_sequence (list): List of DNA bases.
        codon_map (dict): Codon-to-amino-acid mapping.
    Returns:
        list: Protein sequence.
    """
    protein_sequence = []
    for i in range(0, len(dna_sequence) - 2, 3):  # Iterate through triplets
        codon = "".join(dna_sequence[i : i + 3])  # Join triplet bases into a codon
        protein_sequence.append(codon_map.get(codon, "Unknown"))
    return protein_sequence


def find_most_frequent_base(sequence):
    """
    Finds the most frequent base in a DNA sequence.

    Args:
        sequence (list): List of DNA bases.
    Returns:
        str: Most frequent base.
    """
    from collections import Counter

    count = Counter(sequence)
    return count.most_common(1)[0][0]


# Uncomment to run demos
# print(rna_to_dna_conversion(["A", "U", "G", "C", "U"]))
# print(translate_dna_to_protein(["A", "T", "G", "T", "A", "A"], codon_to_amino_acid))
# print(find_most_frequent_base(["A", "T", "G", "C", "A", "A", "G"]))
