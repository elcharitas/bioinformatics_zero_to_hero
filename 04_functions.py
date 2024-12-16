# Module 4: Functions and Modules

"""
Learning Objectives:
- Understand how to write and call functions in Python
- Learn about scope and arguments in functions
- Practice importing and using Python modules

Key Concepts Covered:
1. Writing and Calling Functions
2. Scope and Arguments
3. Importing and Using Python Modules
"""

# 1. Writing and Calling Functions
# Functions allow us to encapsulate code for reuse and organization


def greet_user(name):
    """
    Greets the user by their name

    Args:
        name (str): The name of the user to greet
    Returns:
        None
    """
    print(f"Hello, {name}! Welcome to Module 4.")


# Calling the function
greet_user("Alice")


# Example of a function with a default argument
def greet_with_default(name="User"):
    """
    Greets the user with a default name if no argument is provided

    Args:
        name (str): The name of the user to greet
    Returns:
        None
    """
    print(f"Hello, {name}! How are you today?")


# Calling the function with and without arguments
greet_with_default()
greet_with_default("Bob")


# Example of a function that returns a value
def add_numbers(a, b):
    """
    Adds two numbers and returns the result

    Args:
        a (int or float): First number
        b (int or float): Second number
    Returns:
        int or float: The sum of the two numbers
    """
    return a + b


# Using the function and printing the result
result = add_numbers(5, 7)
print(f"The sum is: {result}")

# 2. Scope and Arguments
# Variables can have global or local scope

global_variable = "I am global!"  # Global variable


def demonstrate_scope(local_argument):
    """
    Demonstrates the concept of variable scope

    Args:
        local_argument (str): A local variable passed to the function
    Returns:
        None
    """
    local_variable = "I am local!"  # Local variable
    print(global_variable)  # Accessing global variable
    print(local_variable)  # Accessing local variable
    print(local_argument)  # Accessing argument


# Calling the function
demonstrate_scope("I am an argument!")


# Example of a function with keyword arguments
def describe_sample(species, count=1, is_pathogenic=False):
    """
    Prints a description of a biological sample

    Args:
        species (str): Name of the species
        count (int, optional): Number of samples. Defaults to 1.
        is_pathogenic (bool, optional): Whether the species is pathogenic. Defaults to False.
    Returns:
        None
    """
    print(f"Species: {species}")
    print(f"Count: {count}")
    print(f"Pathogenic: {is_pathogenic}")


# Calling the function with different arguments
describe_sample("E. coli", 5, True)
describe_sample("S. cerevisiae")

# 3. Importing and Using Python Modules
# Python modules allow us to reuse existing code and functionality

import math  # Importing the math module

# Using a function from the math module
number = 16
square_root = math.sqrt(number)
print(f"The square root of {number} is {square_root}.")

from random import choice, randint  # Importing specific functions from random module

dna_bases = ["A", "T", "C", "G"]
random_base = choice(dna_bases)
print(f"Randomly selected DNA base: {random_base}")

# Example: Simulating a dice roll
dice_roll = randint(1, 6)
print(f"Dice rolled: {dice_roll}")


# Example: Calculating the area of a circle
def circle_area(radius):
    """
    Calculates the area of a circle given its radius

    Args:
        radius (float): The radius of the circle
    Returns:
        float: The area of the circle
    """
    return math.pi * radius**2


# Calculating and printing the area
radius = 3
print(f"The area of a circle with radius {radius} is {circle_area(radius)}")

# Practical Exercise Questions
"""
1. Function Practice:
   Q1: Write a function that calculates the factorial of a given number.
       Example Input: 5
       Expected Output: 120

   Q2: Write a function that takes a list of numbers and returns their average.
       Example Input: [10, 20, 30, 40]
       Expected Output: 25.0

   Q3: Write a function that checks whether a given string is a palindrome.
       Example Input: "radar"
       Expected Output: True

2. Module Practice:
   Q4: Use the random module to simulate the roll of two six-sided dice and print their sum.
   Q5: Use the math module to calculate the circumference of a circle given its radius.
"""


def factorial(n):
    """
    Calculates the factorial of a number

    Args:
        n (int): The number to calculate the factorial for
    Returns:
        int: The factorial of the number
    """
    if n == 0 or n == 1:
        return 1
    else:
        return n * factorial(n - 1)


def calculate_average(numbers):
    """
    Calculates the average of a list of numbers

    Args:
        numbers (list): A list of numeric values
    Returns:
        float: The average of the numbers
    """
    return sum(numbers) / len(numbers)


def is_palindrome(s):
    """
    Checks if a given string is a palindrome

    Args:
        s (str): The string to check
    Returns:
        bool: True if the string is a palindrome, False otherwise
    """
    return s == s[::-1]


# Uncomment the following lines to test the functions
# print(factorial(5))
# print(calculate_average([10, 20, 30, 40]))
# print(is_palindrome("radar"))
