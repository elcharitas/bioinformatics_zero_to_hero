# Module 0.1: Python Basics and Data Types

"""
Learning Objectives:
- Understand fundamental Python data types
- Learn variable assignment and type conversion
- Practice basic data type manipulations

Key Concepts Covered:
1. Variables and Basic Data Types
2. Type Conversion
3. Basic Operators
"""

# 1. Variables and Basic Data Types
# Integer and Float Examples
# Demonstrates how integers and floats work in Python
age = 25  # Integer representing a person’s age
height = 5.9  # Float representing a person’s height in feet
print(f"Age: {age}, Type: {type(age)}")
print(f"Height: {height}, Type: {type(height)}")

# String Manipulation
# Strings are sequences of characters enclosed in quotes
first_name = "Bianca"  # A string representing a first name
last_name = "Rodriguez"  # A string representing a last name
# Concatenating strings to form a full name
full_name = first_name + " " + last_name
print(f"Full Name: {full_name}")

# 2. List Operations
# Lists are ordered, mutable collections of items
# Example: A list of DNA bases
dna_bases = ["A", "T", "C", "G"]  # Initial list of bases
# Adding a new base to the list
# The append() method allows dynamic list modification
dna_bases.append("T")
print(f"DNA Bases: {dna_bases}")

# 3. Dictionary Example
# Dictionaries store key-value pairs
student = {
    "name": "Alex Chen",  # Key: 'name', Value: 'Alex Chen'
    "major": "Bioinformatics",  # Key: 'major', Value: 'Bioinformatics'
    "gpa": 3.8,  # Key: 'gpa', Value: 3.8
}
print(f"Student Info: {student}")

# 4. Type Conversion
# Converting data types using built-in functions
str_number = "42"  # A string representation of a number
# Converting string to integer
int_number = int(str_number)
# Converting string to float
float_number = float(str_number)
print(f"String to Integer: {int_number}")
print(f"String to Float: {float_number}")

# 5. Basic Operators
# Arithmetic Operators
x = 10
y = 3
print(f"Addition: {x} + {y} = {x + y}")
print(f"Subtraction: {x} - {y} = {x - y}")
print(f"Multiplication: {x} * {y} = {x * y}")
print(f"Division: {x} / {y} = {x / y}")
print(f"Modulus: {x} % {y} = {x % y}")

# Comparison Operators
print(f"Is x equal to y? {x == y}")
print(f"Is x not equal to y? {x != y}")
print(f"Is x greater than y? {x > y}")
print(f"Is x less than or equal to y? {x <= y}")

# Logical Operators
is_dna = "A" in dna_bases  # Checks if 'A' is in the list of DNA bases
is_rna = "U" in dna_bases  # Checks if 'U' is in the list of DNA bases
print(f"Is it a DNA sequence? {is_dna}")
print(f"Is it an RNA sequence? {is_rna}")

# Practical Exercise Questions
"""
1. Data Type Conversion Challenge:
   Q1: Write a script that converts a list of string numbers 
       into a list of integers and calculates their sum.
   Example Input: ['10', '20', '30', '40']
   Expected Output: Sum of integers = 100

2. Mixed Data Type Manipulation:
   Q2: Create a dictionary representing a biological sample 
       with at least 3 different data types. Print out each 
       value and its corresponding type.

3. Operator Challenge:
   Q3: Write a script that uses arithmetic operators to calculate 
       the square and cube of a number.
   Q4: Write a script that checks if a given nucleotide is present 
       in a DNA sequence using logical operators.
"""


def type_conversion_challenge(string_numbers):
    """
    Converts list of string numbers to integers and returns sum

    Args:
        string_numbers (list): List of number strings
    Returns:
        int: Sum of converted integers
    """
    # Using list comprehension to convert strings to integers
    int_numbers = [int(num) for num in string_numbers]
    # Returning the sum of integers
    return sum(int_numbers)


def sample_dict_demonstration():
    """
    Creates a dictionary with mixed data types

    Returns:
        dict: Sample biological dictionary
    """
    sample = {
        "species": "E. coli",  # String representing species name
        "colony_count": 1500,  # Integer for colony count
        "is_pathogenic": False,  # Boolean for pathogenicity
        "growth_rate": 0.75,  # Float representing growth rate
    }
    # Iterating through dictionary to print values and their types
    for key, value in sample.items():
        print(f"{key}: {value} (Type: {type(value)})")
    return sample


# Uncomment the following lines to test the functions
# print(type_conversion_challenge(['10', '20', '30', '40']))
# sample_dict_demonstration()
