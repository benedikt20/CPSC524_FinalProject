import subprocess
import os
import csv

# Create the output directory if it doesn't exist
output_dir = "graphs"
os.makedirs(output_dir, exist_ok=True)

# Define the test matrix as (label, N, p)
test_cases = [
    ("A", 100, 0.05),
    ("B", 100, 0.2),
    ("C", 1000, 0.01),
    ("D", 1000, 0.1),
    ("E", 10000, 0.001),
    ("F", 10000, 0.01),
    ("G", 50000, 0.0002),
    ("H", 50000, 0.002)
]

# Path to the compiled ergen binary
ergen_path = "./ergen"

# Prepare the CSV file to log the metrics
csv_filename = os.path.join(output_dir,"metrics.csv")
csv_header = ["Label", "Nodes", "Edges", "Average Degree", "Diameter", "Clustering Coefficient", "File"]

# Check if the file exists to decide whether to write header
file_exists = os.path.isfile(csv_filename)

with open(csv_filename, mode='a', newline='') as csvfile:
    writer = csv.writer(csvfile)
    
    if not file_exists:
        writer.writerow(csv_header)

    for label, n, p in test_cases:
        # Create output filename
        p_str = str(p).replace('.', '')
        filename = f"graph_{label}_n{n}_p{p_str}.txt"
        filepath = os.path.join(output_dir, filename)

        print(f"Generating graph {label}: N={n}, p={p} -> {filepath}")
        
        # Run ergen and capture output
        result = subprocess.run(
            [ergen_path, str(n), str(p), filepath],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Parse the special metrics line
        metrics_line = next((line for line in result.stdout.splitlines() if line.startswith("[METRICS_CSV]")), None)

        if metrics_line:
            parts = metrics_line.split(",")
            if len(parts) == 7:
                _, nodes, edges, avg_deg, diameter, clustering, graph_file = parts
                writer.writerow([label, int(nodes), int(edges), float(avg_deg), int(diameter), float(clustering), graph_file])
            else:
                print(f"Warning: Malformed METRICS_CSV line for {label}")
        else:
            print(f"Warning: No METRICS_CSV found in output for {label}")
