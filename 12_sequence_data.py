# Module 12: Working with Sequence Data

"""
Learning Objectives:
- Learn how to read and write sequence files (FASTA, GenBank).
- Understand the use of `Seq` and `SeqRecord` objects for sequence manipulation.

Key Concepts Covered:
1. Reading and Writing Sequence Files
2. Sequence Manipulation with `Seq` Objects
3. Working with `SeqRecord` Objects
4. Advanced Sequence Parsing and Writing Techniques
"""

# Importing BioPython modules
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# 1. Reading and Writing Sequence Files
# Reading a FASTA file and printing its sequences
def read_fasta(file_path):
    """
    Reads a FASTA file and prints sequence IDs and sequences

    Args:
        file_path (str): Path to the FASTA file
    """
    print(f"Reading FASTA file from: {file_path}")
    for record in SeqIO.parse(file_path, "fasta"):
        print(f"ID: {record.id}")
        print(f"Sequence: {record.seq}")
        print(f"Sequence Length: {len(record.seq)}")


# Writing sequences to a new FASTA file
def write_fasta(file_path, records):
    """
    Writes a list of SeqRecord objects to a FASTA file

    Args:
        file_path (str): Path to the output FASTA file
        records (list): List of SeqRecord objects
    """
    print(f"Writing {len(records)} sequences to FASTA file at: {file_path}")
    SeqIO.write(records, file_path, "fasta")


# Advanced example: Filtering sequences based on length
def filter_fasta_by_length(input_path, output_path, min_length):
    """
    Filters sequences in a FASTA file by length and writes them to a new file

    Args:
        input_path (str): Path to input FASTA file
        output_path (str): Path to output FASTA file
        min_length (int): Minimum length of sequences to retain
    """
    records = [
        record
        for record in SeqIO.parse(input_path, "fasta")
        if len(record.seq) >= min_length
    ]
    write_fasta(output_path, records)


# 2. Sequence Manipulation with `Seq` Objects
# Demonstrating basic sequence manipulations
def seq_operations():
    """
    Demonstrates creating and manipulating sequences with Seq objects
    """
    dna_seq = Seq("ATGCGTACGTT")
    print(f"Original Sequence: {dna_seq}")
    print(f"Length: {len(dna_seq)}")
    print(f"Complement: {dna_seq.complement()}")
    print(f"Reverse Complement: {dna_seq.reverse_complement()}")
    print(f"Transcription to RNA: {dna_seq.transcribe()}")
    print(f"Translation to Protein: {dna_seq.translate(to_stop=True)}")


# 3. Working with `SeqRecord` Objects
# Creating and manipulating a SeqRecord object
def seqrecord_example():
    """
    Demonstrates the use of SeqRecord objects
    """
    # Creating a SeqRecord object
    dna_seq = Seq("ATGCGTACGTT")
    record = SeqRecord(dna_seq, id="seq1", description="Example DNA sequence")
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Sequence: {record.seq}")
    print(f"Sequence Length: {len(record.seq)}")


# Advanced example: Adding features to SeqRecord
def seqrecord_with_features():
    """
    Creates a SeqRecord with annotated features
    """
    dna_seq = Seq("ATGCGTACGTT")
    record = SeqRecord(dna_seq, id="seq2", description="Annotated sequence")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["source"] = "Simulated Data"
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Annotations: {record.annotations}")


# Practical Exercise Questions
"""
1. Reading Sequence Data:
   Q1: Write a function that reads a GenBank file and prints out all the sequence IDs.

2. Writing Sequence Data:
   Q2: Create three SeqRecord objects and write them to a FASTA file. Include IDs and descriptions for each.

3. Sequence Manipulation:
   Q3: Write a script that reads a DNA sequence from a FASTA file and prints its reverse complement.

4. Filtering and Processing:
   Q4: Implement a function to filter sequences from a FASTA file based on a specified length threshold.
"""


def read_genbank_ids(file_path):
    """
    Reads a GenBank file and prints sequence IDs

    Args:
        file_path (str): Path to the GenBank file
    """
    print(f"Reading GenBank file from: {file_path}")
    for record in SeqIO.parse(file_path, "genbank"):
        print(f"ID: {record.id}")


# Uncomment the following lines to test the functions
# seq_operations()
# seqrecord_example()
# seqrecord_with_features()
# read_fasta("example.fasta")
# filter_fasta_by_length("example.fasta", "filtered_output.fasta", 100)
# read_genbank_ids("example.gb")
