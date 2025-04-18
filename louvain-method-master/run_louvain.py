import os
import subprocess
import time
import re

# Define the directory containing the bin files
graphs_dir = './graphs'
main_dir = './' # Main directory where output files will be saved

# Summary file to store metrics
summary_file = os.path.join(main_dir, "louvain_summary_par.txt")

# Get all .bin files in the graphs directory
input_files = [f for f in os.listdir(graphs_dir) if f.endswith('_input_file.bin')]
weight_files = [f for f in os.listdir(graphs_dir) if f.endswith('_weights_file.bin')]

# Output directory for output_community .txt files (results of Lovains algorithm)
output_dir = os.path.join(main_dir, 'output_communities')
os.makedirs(output_dir, exist_ok=True)

# Print the files for debugging
# print("Input files found:", input_files)
# print("Weight files found:", weight_files)

# Create/overwrite summary file with header
with open(summary_file, 'w') as f:
    f.write("Graph Name\tNodes\tLinks\tModularity\tRuntime (s)\n")

# Ensure that each input file has a corresponding weight file
for input_file in input_files:
    # Extract the base name (everything before *input*file.bin)
    base_name = input_file.split('_input_file')[0]
    
    # Find the corresponding weight file
    weight_file = f"{base_name}_weights_file.bin"
    
    if weight_file in weight_files:
        # Define the output file name in the main directory
        output_file = os.path.join(output_dir, f"{base_name}_output_communities.txt")
        
        # Construct the command from the compiled script communityPar.cpp
        command = ['./community_par', os.path.join(graphs_dir, input_file),
                  '-w', os.path.join(graphs_dir, weight_file),
                  '-v', '-o', output_file]
        
        # Start timing
        start_time = time.time()
        
        # Run the command with output piped to a temporary file
        temp_output_file = f"{base_name}_temp_output.txt"
        with open(temp_output_file, 'w') as f_out:
            process = subprocess.run(command, stdout=f_out, stderr=subprocess.STDOUT)
        
        # End timing
        end_time = time.time()
        runtime = end_time - start_time
        
        # Read the output from the temporary file
        with open(temp_output_file, 'r') as f_in:
            output = f_in.read()
        
        # Extract modularity (last numerical value in the output)
        modularity = None
        modularity_match = re.search(r'([\d.]+)\s*$', output)
        if modularity_match:
            modularity = float(modularity_match.group(1))
        
        # Extract number of nodes and links from the first level
        nodes = None
        links = None
        network_info = re.search(r'network size: (\d+) nodes, (\d+) links', output)
        if network_info:
            nodes = network_info.group(1)
            links = network_info.group(2)

        # Print the original output to maintain the same console feedback
        print(output)
        
        # Print summary results
        print(f"Processed: {input_file} and {weight_file} -> {output_file}")
        print(f"Runtime: {runtime:.2f} seconds")
        print(f"Final modularity: {modularity}")
        print(f"Nodes: {nodes}, Links: {links}")
        print("-" * 50)
        
        # Save results to summary file
        with open(summary_file, 'a') as f:
            f.write(f"{base_name}\t{nodes}\t{links}\t{modularity}\t{runtime:.2f}\n")
        
        # Clean up temporary output file
        os.remove(temp_output_file)
    else:
        print(f"Warning: No weight file found for {input_file}")

print(f"Summary of all runs saved to {summary_file}")