# Module 14: Handling Biological Databases

"""
Learning Objectives:
- Learn to access NCBI databases using Bio.Entrez
- Retrieve biological sequence data (e.g., protein or gene sequences)
- Save and parse retrieved sequence data
- Handle potential errors when querying databases

Key Concepts Covered:
1. Introduction to Bio.Entrez
2. Accessing and querying NCBI databases
3. Retrieving sequence data and parsing results
4. Saving and handling sequence data
"""

# 1. Introduction to Bio.Entrez
# Bio.Entrez is a BioPython module that allows access to NCBI databases.
from Bio import Entrez

# Set up an email address to be used for querying NCBI databases
# This is required by NCBI to identify the user and prevent misuse
Entrez.email = "your_email@example.com"  # Replace with your actual email

# 2. Accessing and Querying NCBI Databases
# Example: Searching for a gene (e.g., BRCA1) in the Nucleotide database


def search_gene_in_ncbi(search_term, database="nucleotide", retmax=10):
    """
    Searches for a term in the specified NCBI database.

    Args:
        search_term (str): Search term (e.g., gene name or protein name).
        database (str): The NCBI database to query (default is "nucleotide").
        retmax (int): Maximum number of results to retrieve (default is 10).

    Returns:
        dict: Search results containing record IDs and other metadata.
    """
    try:
        handle = Entrez.esearch(db=database, term=search_term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        return record
    except Exception as e:
        print(f"Error accessing NCBI: {e}")
        return {}


search_term = "BRCA1 Homo sapiens"  # Query for BRCA1 gene in humans
record = search_gene_in_ncbi(search_term)

if record:
    # Print search results
    print(f"Search results for '{search_term}':")
    print(f"Number of records found: {record['Count']}")
    print(f"Record IDs: {record['IdList']}")

# 3. Retrieving Sequence Data
# Retrieve the sequence data for the first record ID in the search results


def fetch_and_save_sequence(
    record_id, database="nucleotide", filename="sequence.fasta"
):
    """
    Fetches sequence data and saves it to a file.

    Args:
        record_id (str): The record ID to fetch data for.
        database (str): The NCBI database to fetch from (default is "nucleotide").
        filename (str): The file to save the retrieved sequence (default is "sequence.fasta").

    Returns:
        str: The retrieved sequence data.
    """
    try:
        handle = Entrez.efetch(
            db=database, id=record_id, rettype="fasta", retmode="text"
        )
        sequence_data = handle.read()
        handle.close()

        # Save to file
        with open(filename, "w") as file:
            file.write(sequence_data)

        print(f"Sequence saved to {filename}")
        return sequence_data
    except Exception as e:
        print(f"Error retrieving sequence: {e}")
        return ""


if record and record["IdList"]:
    record_id = record["IdList"][0]  # Select the first record
    sequence_data = fetch_and_save_sequence(record_id)

    # Print the retrieved sequence
    if sequence_data:
        print(f"\nSequence data for record ID {record_id}:")
        print(sequence_data)
else:
    print("No records found or failed to retrieve data.")

# Practical Exercise Questions
"""
1. NCBI Search and Retrieve:
   Q1: Write a script to search the NCBI Protein database for a specific protein (e.g., hemoglobin) 
       and print the number of records found.

   Q2: Retrieve the FASTA sequence of a protein record using its ID and save it to a file.

2. Parsing Results:
   Q3: Modify the script to extract and print the length of the retrieved nucleotide sequence.

   Q4: Write a script that searches for a gene in a different organism (e.g., TP53 in mice) and prints
       the first three record IDs.

3. Advanced Handling:
   Q5: Implement error handling for scenarios where no records are found or the database query fails.

   Q6: Write a script to fetch multiple sequences by their IDs and save each in a separate file.
"""

# Additional Examples
# Searching for protein sequences in NCBI's protein database
protein_search_term = "hemoglobin Homo sapiens"
protein_record = search_gene_in_ncbi(protein_search_term, database="protein")

if protein_record:
    print(f"\nProtein search results for '{protein_search_term}':")
    print(f"Number of records found: {protein_record['Count']}")
    print(f"Record IDs: {protein_record['IdList']}")

    # Retrieve the first protein sequence and save to a file
    if protein_record["IdList"]:
        protein_id = protein_record["IdList"][0]
        protein_sequence = fetch_and_save_sequence(
            protein_id, database="protein", filename="protein_sequence.fasta"
        )
        print(f"\nProtein sequence data for record ID {protein_id}:")
        print(protein_sequence)
