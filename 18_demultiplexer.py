# Module 18: Demultiplexer CapStone Project

"""Learning Objectives:

- Develop a demultiplexing algorithm for NGS data.
- Implement a Python script for demultiplexing.
- Apply quality control and read trimming techniques.
- Analyze and visualize demultiplexed data.
- Create a comprehensive demultiplexing workflow.
- Perform a hands-on project using real NGS data.
"""

# 1. Introduction to Demultiplexing

# Demultiplexing is the process of separating mixed samples in Next-Generation Sequencing (NGS) data based on unique barcode sequences. This step is essential for downstream analysis and sample identification.

# 2. Project Overview

# In this capstone project, you will develop a demultiplexing algorithm to process NGS data containing mixed samples with unique barcodes. The algorithm will separate reads into individual samples based on barcode matching. You will implement the algorithm in Python and apply quality control and read trimming techniques to ensure data accuracy.
# Check out this sample demultiplexer - https://github.com/jfjlaros/demultiplex/blob/master/demultiplex/demultiplex.py

from Bio import SeqIO


# Function to demultiplex reads based on barcodes
def demultiplex_reads(input_fastq, barcodes, output_files):
    """
    Demultiplexes a FASTQ file based on barcodes.

    Args:
        input_fastq (str): Path to the input FASTQ file.
        barcodes (dict): Dictionary mapping barcodes to sample names.
        output_files (dict): Dictionary mapping sample names to output file paths.

    Returns:
        None
    """
    # Open output files for each sample
    output_handles = {sample: open(path, "w") for sample, path in output_files.items()}

    # Parse the FASTQ file and classify reads
    for record in SeqIO.parse(input_fastq, "fastq"):
        read_sequence = str(record.seq)
        for barcode, sample in barcodes.items():
            if read_sequence.startswith(barcode):
                SeqIO.write(record, output_handles[sample], "fastq")
                break

    # Close all output files
    for handle in output_handles.values():
        handle.close()

