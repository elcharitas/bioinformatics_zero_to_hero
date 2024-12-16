# Module 6: File Handling

"""
Learning Objectives:
- Understand how to read and write files in Python.
- Work with specific bioinformatics file formats like FASTA and CSV.
- Implement error handling in file operations to manage unexpected issues.

Key Concepts Covered:
1. Reading and Writing Files
2. Working with FASTA and CSV Files
3. Error Handling in File Operations
"""

# 1. Reading and Writing Files
# Demonstrates how to read from and write to text files in Python


def read_file(filepath):
    """
    Reads a file and prints its content.

    Args:
        filepath (str): Path to the file to be read.
    """
    try:
        with open(filepath, "r") as file:
            content = file.read()
            print("File Content:")
            print(content)
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")


def write_file(filepath, content):
    """
    Writes content to a file.

    Args:
        filepath (str): Path to the file to be written to.
        content (str): Content to write to the file.
    """
    try:
        with open(filepath, "w") as file:
            file.write(content)
            print(f"Content written to '{filepath}' successfully.")
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")


# 2. Working with FASTA and CSV Files
# FASTA File Example


def read_fasta(filepath):
    """
    Reads a FASTA file and prints each sequence.

    Args:
        filepath (str): Path to the FASTA file.
    """
    try:
        with open(filepath, "r") as file:
            sequence = ""
            for line in file:
                if line.startswith(">"):
                    if sequence:
                        print(f"Sequence: {sequence}")
                    print(f"Header: {line.strip()}")
                    sequence = ""
                else:
                    sequence += line.strip()
            if sequence:
                print(f"Sequence: {sequence}")
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")


# CSV File Example
import csv


def read_csv(filepath):
    """
    Reads a CSV file and prints its rows.

    Args:
        filepath (str): Path to the CSV file.
    """
    try:
        with open(filepath, "r") as file:
            reader = csv.reader(file)
            for row in reader:
                print(row)
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")


# 3. Error Handling in File Operations
# Demonstrates robust file handling with try-except blocks


def robust_file_operations(filepath):
    """
    Example of handling file errors gracefully.

    Args:
        filepath (str): Path to the file for reading or writing.
    """
    try:
        with open(filepath, "r") as file:
            print(file.read())
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
    except IOError:
        print(f"Error: An I/O error occurred with '{filepath}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# Practical Exercise Questions
"""
1. Reading and Writing Files:
   Q1: Write a script to read a text file line by line and count the number of lines.

2. FASTA File Operations:
   Q2: Write a script to extract all sequences from a given FASTA file and save them to a new text file.

3. Error Handling Challenge:
   Q3: Modify the read_csv function to handle empty files gracefully and notify the user if no data is found.
"""

# Uncomment the following lines to test the functions
# write_file('example.txt', 'Hello, Bioinformatics!')
# read_file('example.txt')
# read_fasta('example.fasta')
# read_csv('example.csv')
# robust_file_operations('nonexistent.txt')
