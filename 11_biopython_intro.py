# Module 11: Introduction to BioPython

"""
Learning Objectives:
- Gain an overview of the BioPython library and its applications.
- Learn how to install BioPython.
- Understand the structure and key modules of BioPython.

Key Concepts Covered:
1. What is BioPython?
2. Installing BioPython.
3. Key Modules and Features of BioPython.
"""

# 1. Overview of BioPython
# BioPython is a collection of tools for biological computation.
# It simplifies the process of working with sequence data, file parsing, alignments, and more.
print("What is BioPython?")
print("BioPython is a Python library designed for biological computation.\n")

# 2. Installing BioPython
# Installation can be done via pip.
print("How to Install BioPython?")
print("Run the following command in your terminal to install BioPython:")
print("pip install biopython")

# Example: Importing BioPython to check installation
try:
    import Bio

    print("BioPython installed successfully.")
except ImportError:
    print(
        "BioPython is not installed. Please use the command 'pip install biopython' to install it."
    )

# 3. Key Modules and Features of BioPython
# Overview of some essential BioPython modules
print("Key Modules in BioPython:")
print(
    "1. Bio.Seq: Working with sequences.\n2. Bio.Align: Performing alignments.\n3. Bio.Phylo: Working with phylogenetic trees.\n4. Bio.Entrez: Accessing online biological databases.\n"
)

# Practical Exercise Questions
"""
1. Installation Check:
   Q1: Write a script to check if BioPython is installed on your system.
       If not, display a message guiding the user to install it using pip.

2. Exploring Key Modules:
   Q2: Write a script that imports Bio.Seq, creates a simple DNA sequence, 
       and prints the reverse complement of the sequence.
"""


def check_biopython_installation():
    """
    Checks if BioPython is installed and guides the user if it isn't.

    Returns:
        str: Installation status message.
    """
    try:
        import Bio

        return "BioPython is installed."
    except ImportError:
        return "BioPython is not installed. Run 'pip install biopython' to install it."


def sequence_reverse_complement():
    """
    Creates a DNA sequence and prints its reverse complement using Bio.Seq.

    Returns:
        str: Reverse complement of the DNA sequence.
    """
    from Bio.Seq import Seq

    dna_seq = Seq("ATGC")  # Example DNA sequence
    reverse_complement = dna_seq.reverse_complement()
    return f"Original Sequence: {dna_seq}\nReverse Complement: {reverse_complement}"


# Uncomment the following lines to test the functions
# print(check_biopython_installation())
# print(sequence_reverse_complement())
