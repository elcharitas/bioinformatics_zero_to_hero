# Module 3: Input/Output Operations, Conditional Statements, and Loops

"""
Learning Objectives:
- Master input and output operations
- Learn how to use conditional statements (if-else)
- Understand and implement loops (for and while)
- Apply these concepts to bioinformatics tasks, such as DNA sequence analysis
"""

# 1. Input and Output Operations
# Demonstrates how to take user input and print output
# Example: Getting user input for a DNA sequence
print("Enter a DNA sequence:")
dna_sequence = input()  # User inputs a DNA sequence
print(f"You entered: {dna_sequence}")

# Writing and reading from a file
# Writing the DNA sequence to a file
with open("dna_sequence.txt", "w") as file:
    file.write(dna_sequence)

# Reading the DNA sequence back from the file
with open("dna_sequence.txt", "r") as file:
    saved_sequence = file.read()
    print(f"Saved DNA sequence: {saved_sequence}")

# 2. Conditional Statements (if-else)
# Checking the composition of a DNA sequence
if "A" in dna_sequence:
    print("The sequence contains Adenine (A).")
else:
    print("The sequence does not contain Adenine (A).")

# Identifying if a sequence contains invalid characters
valid_bases = {"A", "T", "C", "G"}
invalid_bases = set(dna_sequence) - valid_bases
if invalid_bases:
    print(f"Invalid bases found: {invalid_bases}")
else:
    print("The DNA sequence is valid.")

# 3. Loops (for and while)
# Using a for loop to count each base in the sequence
base_counts = {"A": 0, "T": 0, "C": 0, "G": 0}
for base in dna_sequence:
    if base in base_counts:
        base_counts[base] += 1
print(f"Base counts: {base_counts}")

# Using a while loop to repeat a task
# Example: Validate user input until a valid DNA sequence is entered
while True:
    print("Enter a valid DNA sequence (only A, T, C, G):")
    dna_sequence = input()
    if set(dna_sequence).issubset(valid_bases):
        print("Valid DNA sequence entered.")
        break
    else:
        print("Invalid sequence. Try again.")

# 4. Bioinformatics Use Case: DNA Reverse Complement
# Generating the reverse complement of a DNA sequence
complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
reverse_complement = "".join(complement[base] for base in reversed(dna_sequence))
print(f"Reverse complement: {reverse_complement}")

# Practical Exercise Questions
"""
1. Input and Output Challenge:
   Q1: Write a script to read a protein sequence from the user and save it to a file. 
       Then, read the file back and print the sequence.

2. Conditional and Loop Challenges:
   Q2: Create a script to calculate the GC content (percentage of G and C bases) 
       of a DNA sequence provided by the user.

3. DNA Reverse Complement Challenge:
   Q3: Write a function that takes a DNA sequence as input and returns its reverse complement.
"""


def gc_content(sequence):
    """
    Calculates the GC content of a DNA sequence.

    Args:
        sequence (str): DNA sequence
    Returns:
        float: GC content percentage
    """
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def reverse_complement(sequence):
    """
    Generates the reverse complement of a DNA sequence.

    Args:
        sequence (str): DNA sequence
    Returns:
        str: Reverse complement sequence
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in reversed(sequence))


# Uncomment to test the functions
# print(gc_content("ATGCGCTA"))
# print(reverse_complement("ATGCGCTA"))
