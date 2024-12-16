# Module 15: Phylogenetics with BioPython

"""
Learning Objectives:
- Learn to parse phylogenetic trees using Bio.Phylo.
- Understand how to visualize and manipulate phylogenetic trees.

Key Concepts Covered:
1. Parsing Phylogenetic Trees
2. Tree Visualization
3. Tree Manipulation
"""

# Import necessary modules
from Bio import Phylo
import matplotlib.pyplot as plt

# 1. Parsing Phylogenetic Trees
# Parsing a Newick-formatted phylogenetic tree file
# Example: Load a tree from a sample file
example_tree_file = "example_tree.nwk"  # Replace with your Newick file path
# Using Phylo.read() to parse the tree
tree = Phylo.read(example_tree_file, "newick")
print("Parsed Tree:")
Phylo.draw_ascii(tree)  # Display the tree in ASCII format


# Parsing multiple trees from a file
def parse_multiple_trees(file_path):
    """
    Parses multiple trees from a Newick file.

    Args:
        file_path (str): Path to the Newick file containing multiple trees.
    Returns:
        list: List of phylogenetic tree objects.
    """
    trees = list(Phylo.parse(file_path, "newick"))
    print(f"Parsed {len(trees)} trees from {file_path}")
    return trees


# Uncomment the following line to parse multiple trees
# trees = parse_multiple_trees("multiple_trees.nwk")


# 2. Tree Visualization
# Graphical visualization of the tree
def visualize_tree(tree, title="Phylogenetic Tree"):
    """
    Visualizes a phylogenetic tree using matplotlib.

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree object.
        title (str): Title of the tree visualization.
    """
    Phylo.draw(tree)
    plt.title(title)
    plt.show()


# Uncomment to visualize the parsed tree
# visualize_tree(tree, title="Example Phylogenetic Tree")


# Advanced Visualization with annotations
def advanced_visualize_tree(tree):
    """
    Visualizes a tree with advanced annotations, including branch labels.

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree object.
    """
    Phylo.draw(tree, do_show=False)
    plt.title("Annotated Phylogenetic Tree")
    plt.show()


# 3. Tree Manipulation
# Example: Accessing and modifying tree nodes
print("Clade names in the tree:")
for clade in tree.find_clades():
    print(clade.name)


# Pruning a clade
def prune_clade(tree, clade_name):
    """
    Removes a clade from the phylogenetic tree.

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree object.
        clade_name (str): Name of the clade to prune.
    """
    clade_to_prune = tree.find_any(name=clade_name)
    if clade_to_prune:
        tree.prune(clade_to_prune)
        print(f"Pruned clade: {clade_name}")
    else:
        print(f"Clade {clade_name} not found.")


# Uncomment to prune a clade
# prune_clade(tree, "CladeA")


# Counting leaf nodes
def count_leaf_nodes(tree):
    """
    Counts the number of leaf nodes in the given phylogenetic tree.

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree object
    Returns:
        int: Count of leaf nodes
    """
    leaf_count = 0
    for clade in tree.find_clades():
        if not clade.clades:  # A clade with no descendants is a leaf
            leaf_count += 1
    return leaf_count


# Uncomment to count leaf nodes
# print("Number of leaf nodes:", count_leaf_nodes(tree))


# Saving the modified tree to a file
def save_tree(tree, file_path):
    """
    Saves the modified tree to a Newick file.

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree object.
        file_path (str): Path to save the tree.
    """
    Phylo.write(tree, file_path, "newick")
    print(f"Tree saved to {file_path}")


# Uncomment to save the modified tree
# save_tree(tree, "modified_tree.nwk")

# Practical Exercise Questions
"""
1. Parsing Trees:
   Q1: Write a script to parse and display a phylogenetic tree 
       in ASCII format from a given Newick file.

2. Tree Visualization:
   Q2: Visualize a phylogenetic tree using Bio.Phylo and matplotlib. 
       Add a title to the visualization.

3. Tree Manipulation:
   Q3: Write a script to iterate through all clades in a tree and 
       count the number of clades with no descendants (leaf nodes).

4. Advanced Manipulation:
   Q4: Write a script to rename a specific clade and save 
       the updated tree to a new file.
"""
