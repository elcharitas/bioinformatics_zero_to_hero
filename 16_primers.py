# Module 16: Primer and Barcode Design with Biopython

"""
Learning Objectives:
- Understand the basics of primer and barcode design.
- Explore PCR principles and primer design considerations.
- Utilize Primer3Plus for automated primer design.
- Implement a simple primer design algorithm in Python.
- Learn the principles of DNA barcoding and design effective barcode sequences.

Key Concepts Covered:
1. Primer Design Basics
2. PCR Principles
3. Primer Design Considerations
4. Barcode Design Principles
5. Hands-on Project: Multiplexed PCR Primer and Barcode Design
"""

# 1. Primer Design Basics
# Primers are short DNA sequences that provide a starting point for DNA synthesis.
# Effective primer design is crucial for successful PCR amplification.

# Example: A simple primer sequence
forward_primer = "ATGCGTACGTCAGT"  # Forward primer sequence
reverse_primer = "ACTGACGTACGCAT"  # Reverse primer sequence
print(f"Forward Primer: {forward_primer}")
print(f"Reverse Primer: {reverse_primer}")

# 2. Primer Design Considerations
# Factors include length, GC content, and melting temperature (Tm).
# Example: Calculating GC content


def calculate_gc_content(sequence):
    """
    Calculates the GC content of a given DNA sequence.

    Args:
        sequence (str): DNA sequence
    Returns:
        float: GC content percentage
    """
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


# Example usage
primer_gc_content = calculate_gc_content(forward_primer)
print(f"GC Content of Forward Primer: {primer_gc_content:.2f}%")

# 3. Barcode Design Principles
# DNA barcoding uses unique DNA sequences to identify samples.
# Effective barcodes are unique and avoid repetitive regions.

# Example: Generating a simple barcode
barcode = "ACTGACTGACTG"
print(f"Barcode: {barcode}")

# 4. Implementing a Simple Primer Design Algorithm in Python
# Design primers based on GC content and sequence length


def design_primer(sequence, target_gc_content, length):
    """
    Designs a primer with the specified GC content and length.

    Args:
        sequence (str): DNA sequence
        target_gc_content (float): Desired GC content
        length (int): Desired primer length
    Returns:
        str: Designed primer
    """
    for i in range(len(sequence) - length + 1):
        candidate = sequence[i : i + length]
        if abs(calculate_gc_content(candidate) - target_gc_content) <= 5:  # Tolerance
            return candidate
    return "No suitable primer found."


# Example usage
sequence = "ATGCGTACGTCAGTACGCGTACGTCAGTACG"
primer = design_primer(sequence, target_gc_content=50, length=15)
print(f"Designed Primer: {primer}")

# 5. Hands-on Project: Multiplexed PCR Experiment
# Design primers and barcodes for a PCR experiment

"""
Project Steps:
1. Choose target DNA sequences for amplification.
2. Use the provided `design_primer` function to create primers for each target.
3. Generate unique barcodes for each sample.
4. Combine primers and barcodes for multiplexed PCR.
"""


def multiplex_pcr_experiment(target_sequences, barcodes):
    """
    Designs primers and assigns barcodes for multiplexed PCR.

    Args:
        target_sequences (list): List of target DNA sequences
        barcodes (list): List of unique barcodes
    Returns:
        dict: Mapping of target sequences to primers and barcodes
    """
    experiment_plan = {}
    for i, sequence in enumerate(target_sequences):
        primer = design_primer(sequence, target_gc_content=50, length=15)
        experiment_plan[sequence] = {
            "primer": primer,
            "barcode": barcodes[i] if i < len(barcodes) else "No barcode assigned",
        }
    return experiment_plan


# Example usage
targets = ["ATGCGTACGTCAGTACG", "GTCAGTACGCGTACGTT"]
bc = ["ACTGACTG", "GTCAGTCA"]
plan = multiplex_pcr_experiment(targets, bc)
for target, details in plan.items():
    print(
        f"Target: {target}\nPrimer: {details['primer']}\nBarcode: {details['barcode']}\n"
    )
