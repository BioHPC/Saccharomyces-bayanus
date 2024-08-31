#!/usr/bin/env python3

import os
import subprocess

# Set input and output paths
input_dir = "/mnt/b16/06_yeast_ONT/19_orthofinder/cds_for_seven_species/cds_orthofinder/Results_Jun30/Orthogroup_Sequences"
output_dir = "/mnt/b16/06_yeast_ONT/19_orthofinder/cds_for_seven_species/cds_orthofinder/Results_Jun30/Orthogroups/mafft_all"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Iterate over all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".fa"):
        input_file = os.path.join(input_dir, filename)
        output_file = os.path.join(output_dir, f"{filename}.mafft.out")
        
        # Run the MAFFT command
        mafft_command = ["mafft", "--auto", input_file]
        
        with open(output_file, "w") as out_f:
            subprocess.run(mafft_command, stdout=out_f)
        
        print(f"Processed {input_file} -> {output_file}")

print("All files processed.")
