# Module 8: Data Visualization in Python

"""
Learning Objectives:
- Understand the basics of data visualization in Python using matplotlib.
- Create plots to represent bioinformatics data.
- Generate a GC content graph as an example of bioinformatics visualization.

Key Concepts Covered:
1. Introduction to matplotlib
2. Line plots and bar plots
3. GC Content visualization
"""

# 1. Introduction to matplotlib
# Importing the matplotlib library
import matplotlib.pyplot as plt

# Creating a simple line plot to demonstrate matplotlib
# Data for plotting
days = [1, 2, 3, 4, 5]
temperatures = [20, 21, 19, 22, 20]

# Plotting the data
plt.plot(days, temperatures, marker="o")
plt.title("Temperature Over Days")
plt.xlabel("Days")
plt.ylabel("Temperature (°C)")
plt.grid(True)
plt.show()

# 2. Line plots and bar plots
# Demonstrating a bar plot
# Sample nucleotide counts
data = {"A": 30, "T": 25, "G": 20, "C": 15}
nucleotides = list(data.keys())
counts = list(data.values())

# Creating a bar plot
plt.bar(nucleotides, counts, color=["blue", "green", "red", "purple"])
plt.title("Nucleotide Counts")
plt.xlabel("Nucleotides")
plt.ylabel("Counts")
plt.show()

# 3. GC Content Visualization
# Function to calculate and plot GC content


def plot_gc_content(sequence):
    """
    Calculates and plots the GC content of a DNA sequence.

    Args:
        sequence (str): A DNA sequence.
    Returns:
        None
    """
    # Counting G and C nucleotides
    gc_count = sequence.count("G") + sequence.count("C")
    total_length = len(sequence)
    gc_content = (gc_count / total_length) * 100

    # Plotting GC content as a pie chart
    labels = ["GC Content", "Other Bases"]
    sizes = [gc_content, 100 - gc_content]
    colors = ["cyan", "lightgrey"]
    explode = (0.1, 0)  # Highlighting GC content

    plt.pie(
        sizes,
        explode=explode,
        labels=labels,
        colors=colors,
        autopct="%1.1f%%",
        startangle=140,
    )
    plt.title("GC Content of the Sequence")
    plt.axis("equal")  # Equal aspect ratio ensures the pie chart is circular
    plt.show()


# Example Sequence for GC Content Visualization
example_sequence = "ATGCGTACGTAGCTAGCTAGCTAGCGTAGCTAGCGTACGTA"
plot_gc_content(example_sequence)

# Practical Exercise Questions
"""
1. Basic Visualization Challenge:
   Q1: Create a line plot representing the growth of a bacterial colony over 7 days.
       Example Data: [100, 200, 400, 800, 1600, 3200, 6400]

2. GC Content Visualization:
   Q2: Write a script that calculates the GC content for a list of sequences and generates a bar plot 
       comparing the GC content of each sequence.
"""
