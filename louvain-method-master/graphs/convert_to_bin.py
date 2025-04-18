# File: convert_to_bin.py
# -- Converts weighted graphs from .txt files to input.bin and weights.bin files for community.cpp (lovains algorithm)
# -----------------------------------------------------------------------------
# Reads .txt file in the directory with the line format: i j w
# where i and j are nodes with associated weight w
#
# Outputs input.bin and weights.bin files that will be accessed by community.cpp (lovains method script)
#
# Uses `convert` compiled file from main_convert.cpp script (that is compiled with `make` command)
#
# This program must not be distributed without agreement of the above mentionned authors.
# -----------------------------------------------------------------------------
# Author   : Benedikt Farag
# Email    : benedikt.farag@yale.edu
# Location : New Haven CT
# Time	    : April 2025
# -----------------------------------------------------------------------------
# To execute  : python3 convert_to_bin.py 

import os
import subprocess

# Get all .txt files in the current directory
txt_files = [f for f in os.listdir() if f.endswith('.txt')]

# Loop through each .txt file and run the convert command
for txt_file in txt_files:
    # Define the output file names based on the input file
    base_name = os.path.splitext(txt_file)[0]
    input_file = txt_file
    output_file_bin = f"{base_name}_input_file.bin"
    weights_file_bin = f"{base_name}_weights_file.bin"
    
    # Construct the command
    command = ['./convert', '-i', input_file, '-o', output_file_bin, '-w', weights_file_bin]
    
    # Run the command
    subprocess.run(command, check=True)
    print(f"Processed: {input_file} -> {output_file_bin}, {weights_file_bin}")
