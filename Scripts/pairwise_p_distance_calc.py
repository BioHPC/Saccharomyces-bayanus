import os
from Bio import SeqIO

def calculate_p_distance(seq1, seq2):
    differences = 0
    valid_positions = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != '-' and c2 != '-':
            valid_positions += 1
            if c1 != c2:
                differences += 1
    return differences / valid_positions if valid_positions > 0 else 0

# Path settings
input_dir = './mafft_all'
output_path = './output_p_distance.csv'

# Prepare the output file
with open(output_path, 'w') as output_file:
    output_file.write("Filename,Sbay_Gene_Pair,P_Distance\n")

    # Iterate over all files in the mafft_all folder
    for filename in os.listdir(input_dir):
        if filename.endswith("mafft.out"):
            file_path = os.path.join(input_dir, filename)
            sbay_seqs = []
            non_sbay_seqs = []

            # Read sequence files
            for record in SeqIO.parse(file_path, "fasta"):
                if record.id.startswith("Sbay"):
                    sbay_seqs.append(record)
                else:
                    non_sbay_seqs.append(record)

            # Compare each pair of Sbay and non-Sbay genes
            for sbay in sbay_seqs:
                for non_sbay in non_sbay_seqs:
                    p_distance = calculate_p_distance(sbay.seq, non_sbay.seq)
                    pair_name = f"{sbay.id}_vs_{non_sbay.id}"
                    output_file.write(f"{filename},{pair_name},{p_distance:.4f}\n")
                    output_file.flush()
            print(f"Processed file: {filename}")  # Print processing progress

print("Calculation complete, results saved in: output_p_distance.csv")
