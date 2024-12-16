"""
Biopython Fundamentals Module
==============================
Learning Objectives:
- Understand Biopython's core data structures
- Manipulate biological sequences
- Parse sequence files
- Perform basic sequence analysis

Key Libraries to Import:
- Bio.Seq: Sequence manipulation
- Bio.SeqIO: Sequence file parsing
- Bio.SeqRecord: Sequence record handling
"""

from Bio.Seq import Seq
from Bio import SeqIO
import random


class BiopythonSequenceAnalyzer:
    """
    Comprehensive class for demonstrating Biopython sequence operations
    """

    def __init__(self, sequence):
        """
        Initialize sequence with Biopython Seq object

        Args:
            sequence (str): DNA/RNA sequence
        """
        self.seq = Seq(sequence.upper())

    def transcribe(self):
        """
        Transcribe DNA to RNA

        Returns:
            Seq: Transcribed RNA sequence
        """
        return self.seq.transcribe()

    def translate(self):
        """
        Translate DNA to protein sequence

        Returns:
            Seq: Translated protein sequence
        """
        return self.seq.translate()

    def reverse_complement(self):
        """
        Generate reverse complement of DNA sequence

        Returns:
            Seq: Reverse complement sequence
        """
        return self.seq.reverse_complement()

    def gc_content(self):
        """
        Calculate GC content of sequence

        Returns:
            float: GC content percentage
        """
        gc_count = self.seq.count("G") + self.seq.count("C")
        return (gc_count / len(self.seq)) * 100


def sequence_file_parsing_demo():
    """
    Demonstrate Biopython file parsing capabilities

    Practice Questions:
    1. File Parsing Challenge:
       Q1: Write a function to parse a GenBank file and extract
           specific gene sequences and their annotations.

    2. Sequence Filtering:
       Q2: Create a script that filters sequences based on length,
           GC content, and other molecular characteristics.
    """

    # Simulated GenBank file parsing
    def parse_genbank_file(filename):
        """
        Parse GenBank file and extract sequence information

        Args:
            filename (str): Path to GenBank file
        Returns:
            list: Extracted sequence records
        """
        try:
            # In real scenario, replace with actual file parsing
            records = list(SeqIO.parse(filename, "genbank"))

            # Filter and process records
            processed_records = []
            for record in records:
                seq_analyzer = BiopythonSequenceAnalyzer(str(record.seq))
                processed_records.append(
                    {
                        "id": record.id,
                        "description": record.description,
                        "length": len(record.seq),
                        "gc_content": seq_analyzer.gc_content(),
                        "full_sequence": str(record.seq),
                    }
                )

            return processed_records

        except FileNotFoundError:
            print(f"File not found: {filename}")
            return []


def molecular_sequence_generator():
    """
    Generate random molecular sequences for analysis

    Demonstrates sequence generation and manipulation

    Practice Questions:
    1. Sequence Generation Challenge:
       Q1: Extend this function to generate sequences with
           specific GC content or codon bias.

    2. Mutation Simulation:
       Q2: Create a method to introduce random mutations
           into generated sequences.
    """
    bases = ["A", "T", "C", "G"]

    def generate_dna_sequence(length=100):
        """
        Generate random DNA sequence

        Args:
            length (int): Desired sequence length
        Returns:
            Seq: Randomly generated DNA sequence
        """
        sequence = "".join(random.choice(bases) for _ in range(length))
        return Seq(sequence)

    def mutate_sequence(sequence, mutation_rate=0.01):
        """
        Introduce random mutations to a sequence

        Args:
            sequence (Seq): Original sequence
            mutation_rate (float): Probability of mutation
        Returns:
            Seq: Mutated sequence
        """
        mutated_bases = []
        for base in sequence:
            if random.random() < mutation_rate:
                # Replace with random base
                mutated_bases.append(random.choice(bases))
            else:
                mutated_bases.append(base)

        return Seq("".join(mutated_bases))

    # Generate and demonstrate sequences
    original_seq = generate_dna_sequence(200)
    mutated_seq = mutate_sequence(original_seq)

    print("Original Sequence Length:", len(original_seq))
    print("Mutated Sequence Length:", len(mutated_seq))
    print(
        "Mutation Differences:",
        sum(1 for a, b in zip(original_seq, mutated_seq) if a != b),
    )


def advanced_sequence_analysis():
    """
    Advanced sequence analysis techniques

    Demonstrates complex sequence manipulations
    """
    # Example DNA sequence
    dna_seq = Seq("ATGGCCATGGCGCCCAGAACTGAGATGG")

    # Multiple sequence analyses
    analyses = {
        "Original Sequence": dna_seq,
        "Transcribed RNA": dna_seq.transcribe(),
        "Translated Protein": dna_seq.translate(),
        "Reverse Complement": dna_seq.reverse_complement(),
        "GC Content": f"{(dna_seq.count('G') + dna_seq.count('C')) / len(dna_seq) * 100:.2f}%",
    }

    for name, result in analyses.items():
        print(f"{name}: {result}")


# Demonstration functions (commented out for manual execution)
# sequence_file_parsing_demo()
# molecular_sequence_generator()
# advanced_sequence_analysis()


def main():
    """
    Main execution function demonstrating Biopython capabilities
    """
    # Create a sequence analyzer
    analyzer = BiopythonSequenceAnalyzer("ATGGCCATGGCGCCCAGAACTGAGATGG")

    print("Sequence Analysis Results:")
    print(f"Original Sequence: {analyzer.seq}")
    print(f"Transcription: {analyzer.transcribe()}")
    print(f"Translation: {analyzer.translate()}")
    print(f"Reverse Complement: {analyzer.reverse_complement()}")
    print(f"GC Content: {analyzer.gc_content():.2f}%")


if __name__ == "__main__":
    main()
