# Module 1.7: Introduction to Libraries

"""
Learning Objectives:
- Gain an overview of commonly used Python libraries in bioinformatics (NumPy, pandas, matplotlib).
- Learn how to install and import libraries using pip.

Key Concepts Covered:
1. Overview of Libraries
2. Installing Libraries
3. Basic Usage Examples
"""

# 1. Overview of Libraries
# NumPy: Library for numerical computations
# pandas: Library for data manipulation and analysis
# matplotlib: Library for data visualization

print("Libraries Overview:")
print("- NumPy: Numerical computations like arrays and linear algebra.")
print("- pandas: Data manipulation with DataFrames.")
print("- matplotlib: Plotting and visualization.")

# 2. Installing Libraries
# Ensure you have pip installed (Python's package installer).
# Example installation commands (uncomment to run in your environment):
# !pip install numpy
# !pip install pandas
# !pip install matplotlib

print("\nTo install libraries, use commands like:")
print("pip install numpy pandas matplotlib")

# 3. Basic Usage Examples
# Example 1: Using NumPy
import numpy as np

# Creating a NumPy array
array = np.array([1, 2, 3, 4, 5])
print("\nNumPy Array:", array)
print("Array Mean:", np.mean(array))

# Example 2: Using pandas
import pandas as pd

# Creating a DataFrame
data = {"Gene": ["BRCA1", "TP53", "MYC"], "Expression": [5.4, 3.2, 7.8]}
df = pd.DataFrame(data)
print("\nPandas DataFrame:")
print(df)

# Example 3: Using matplotlib
import matplotlib.pyplot as plt

# Creating a simple bar plot
genes = data["Gene"]
expression_levels = data["Expression"]
plt.bar(genes, expression_levels)
plt.title("Gene Expression Levels")
plt.xlabel("Gene")
plt.ylabel("Expression")
plt.show()

# Practical Exercise Questions
"""
1. Library Installation Challenge:
   Q1: Write a script to check if NumPy, pandas, and matplotlib are installed. 
       If not, print a message indicating which library is missing.

2. Basic Usage Challenge:
   Q2: Create a pandas DataFrame with gene names and GC content values.
       Plot the GC content values using matplotlib.
"""


def check_library_installation():
    """
    Checks if NumPy, pandas, and matplotlib are installed.

    Returns:
        list: Missing libraries
    """
    required_libraries = ["numpy", "pandas", "matplotlib"]
    missing_libraries = []
    for lib in required_libraries:
        try:
            __import__(lib)
        except ImportError:
            missing_libraries.append(lib)
    return missing_libraries


def plot_gc_content():
    """
    Creates a DataFrame with gene names and GC content, then plots the data.
    """
    # Example DataFrame
    data = {"Gene": ["Gene1", "Gene2", "Gene3"], "GC_Content": [45, 50, 55]}
    df = pd.DataFrame(data)
    print("\nGC Content DataFrame:")
    print(df)

    # Plotting GC Content
    plt.bar(df["Gene"], df["GC_Content"], color="blue")
    plt.title("GC Content per Gene")
    plt.xlabel("Gene")
    plt.ylabel("GC Content (%)")
    plt.show()


# Uncomment the following lines to test the functions
# missing_libraries = check_library_installation()
# if missing_libraries:
#     print("Missing Libraries:", missing_libraries)
# else:
#     print("All required libraries are installed.")

# plot_gc_content()
