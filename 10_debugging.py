# Module 10: Debugging and Optimization

"""
Learning Objectives:
- Understand debugging techniques and tools to identify and fix errors.
- Learn methods to optimize Python code for performance improvements.

Key Concepts Covered:
1. Debugging Techniques
2. Using Debugging Tools
3. Code Optimization Strategies
"""


# 1. Debugging Techniques
# Debugging helps identify and fix errors in code.
def debugging_example():
    """
    Example demonstrating debugging techniques
    """
    # Introduce a bug: Division by zero
    try:
        result = 10 / 0  # Intentional bug
    except ZeroDivisionError as e:
        print(f"Error encountered: {e}")
    finally:
        print("Debugging Example Complete")


def logging_example():
    """
    Using logging for debugging and tracing code
    """
    import logging

    logging.basicConfig(level=logging.DEBUG)

    # Debugging a function with logging
    def add_numbers(a, b):
        logging.debug(f"Adding {a} and {b}")
        return a + b

    result = add_numbers(5, 10)
    logging.info(f"Result: {result}")


# 2. Using Debugging Tools
# Demonstrating the use of pdb (Python Debugger)
def pdb_debug_example():
    """
    Example demonstrating the use of pdb for interactive debugging
    """
    import pdb

    x = 5
    y = 10
    pdb.set_trace()  # Start debugging session here
    z = x + y
    print(f"Sum: {z}")


# 3. Code Optimization Strategies
# Demonstrating code optimization techniques
import time


def inefficient_function():
    """
    Example of an inefficient function
    """
    result = 0
    for i in range(10**6):
        result += i
    return result


def optimized_function():
    """
    Optimized version using a built-in function
    """
    return sum(range(10**6))


def measure_execution_time():
    """
    Measure execution time of functions
    """
    start = time.time()
    inefficient_result = inefficient_function()
    end = time.time()
    print(f"Inefficient function took: {end - start:.2f} seconds")

    start = time.time()
    optimized_result = optimized_function()
    end = time.time()
    print(f"Optimized function took: {end - start:.2f} seconds")


# Practical Exercise Questions
"""
1. Debugging Challenge:
   Q1: Write a script that intentionally raises a ValueError. Use try-except-finally to handle it gracefully.
   Example Output: "Error encountered: invalid literal for int() with base 10"

2. Optimization Challenge:
   Q2: Write a script to optimize a function that calculates the factorial of a number using both recursion and iteration. Compare their execution times.
"""


def value_error_debugging():
    """
    Handles ValueError with debugging example
    """
    try:
        number = int("abc")  # Intentional ValueError
    except ValueError as e:
        print(f"Error encountered: {e}")
    finally:
        print("ValueError Debugging Complete")


def factorial_comparison(n):
    """
    Compare recursive and iterative factorial calculations

    Args:
        n (int): Number for factorial calculation
    """

    # Recursive implementation
    def recursive_factorial(num):
        return 1 if num == 0 else num * recursive_factorial(num - 1)

    # Iterative implementation
    def iterative_factorial(num):
        result = 1
        for i in range(1, num + 1):
            result *= i
        return result

    # Measure execution times
    start = time.time()
    recursive_result = recursive_factorial(n)
    end = time.time()
    print(f"Recursive factorial took: {end - start:.6f} seconds")

    start = time.time()
    iterative_result = iterative_factorial(n)
    end = time.time()
    print(f"Iterative factorial took: {end - start:.6f} seconds")

    # Print results
    print(f"Recursive Result: {recursive_result}, Iterative Result: {iterative_result}")


# Uncomment the following lines to test the functions
# debugging_example()
# logging_example()
# pdb_debug_example()
# measure_execution_time()
# value_error_debugging()
# factorial_comparison(10)
