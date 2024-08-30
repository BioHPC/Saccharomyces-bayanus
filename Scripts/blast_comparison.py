import os
import subprocess
import tempfile
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def local_blast_analysis(hybrid_assembly_path, parent_genomes_paths, output_path, chunk_size, score_type, score_threshold):
    # header
    with open(output_path, 'w') as output_file:
        headers = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        output_file.write("\t".join(headers) + "\n")

    with open(hybrid_assembly_path, 'r') as hybrid_assembly_file:
        for record in SeqIO.parse(hybrid_assembly_file, "fasta"):
            sequence = str(record.seq)
            seq_id = record.id 
            current_chunk_sequence = ''
            chunk_start = 1

            for i in range(0, len(sequence), chunk_size):
                chunk_end = chunk_start + chunk_size - 1
                chunk_sequence = sequence[i:i+chunk_size]
                chunk_id = f"{seq_id}_{chunk_start}-{chunk_end}"  

                print("\nChunk: ", chunk_start, " : ", chunk_end)
                perform_local_blast(chunk_sequence, chunk_id, parent_genomes_paths, output_path, score_threshold, score_type)

                chunk_start = chunk_end + 1

    # end (partial) chunk
    if current_chunk_sequence:
        chunk_end = chunk_start + len(current_chunk_sequence) - 1
        perform_local_blast(chunk_to_process, parent_genomes_paths, output_path, score_threshold, args.score_type)

def perform_local_blast(chunk_sequence, chunk_id, parent_genomes_paths, output_path, score_threshold, score_type):
    chunk_record = SeqRecord(Seq(chunk_sequence), id=chunk_id, description="")
    best_scores = []
    best_hits = []
    # Can't pass directly -- throws error due to length
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_seq_file:
        SeqIO.write(chunk_record, temp_seq_file, "fasta")
        temp_seq_file_path = temp_seq_file.name

    for parent_genome_path in parent_genomes_paths:
        try:
            blastn_output = subprocess.check_output(['blastn', '-query', temp_seq_file_path, '-db', parent_genome_path, '-outfmt', '6'])
            decoded_blastn_output = blastn_output.decode('utf-8')
            hits = blastn_output.decode('utf-8').strip().split('\n')
            print("BLAST output for: ", parent_genome_path)
            print(hits[0])
            if hits and hits != ['']:
                best_hit = hits[0].split('\t')
                #if scoretype == 'bitscore':
                score = float(best_hit[11])  # add escore for comparison type (?) 
                #if scoretype == 'escore':
                #score=flat(best_hit[10]):
                best_scores.append(score)
                best_hits.append(best_hit)
            else:
                print("no hits")
                best_scores.append(-1)  # no hits (shouldn't happen)
                best_hits.append(None)

        except subprocess.CalledProcessError as e:
            print("Error: ", e.output)
            best_scores.append(-1)
            best_hits.append(None)

    os.remove(temp_seq_file_path)

    for hit in best_hits:
        if hit:
            #if scoretype=='bitscore':
            score = float(hit[11])  # implement escore (column 10)
            #if scoretype=='escore':
            #score=float(hit[10])
            best_scores.append(score)
        else:
            best_scores.append(-1)

    if best_scores:
        sorted_scores_with_parents = sorted(zip(best_scores, parent_genomes_paths, best_hits), reverse=True)

        if len(sorted_scores_with_parents) > 1 and sorted_scores_with_parents[0][0] - sorted_scores_with_parents[1][0] >= score_threshold:
            inferred_parent = sorted_scores_with_parents[0][1]
            best_hit = sorted_scores_with_parents[0][2]
        else:
            inferred_parent = "Undetermined"
            best_hit = None

        # output
        print("Inferred parent: ",inferred_parent)
        with open(output_path, 'a') as output_file:
            output_line = '\t'.join(best_hit) if best_hit else "No clear hit"
            output_file.write(inferred_parent + "\t" + output_line + "\n")
    else:
        print("No hits found for either parent.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BLAST, output inferred parent based")
    
    parser.add_argument("query_sequence", help="query sequence")
    parser.add_argument("target_sequences", nargs='+', help="target (parental) sequences")
    parser.add_argument("output_path", help="output file")
    parser.add_argument("--chunk_size", type=int, default=5000, help="Chunk size, default is 5k")
    parser.add_argument("--score_type", choices=['bitscore', 'escore'], type=str, default='bitscore', help="'bitscore' or 'escore' (only bitscore currently), bitscore is default")
    parser.add_argument("--score_threshold", type=float, default=100, help="threshold for bitscore/escore (default 100)")

    args = parser.parse_args()

    local_blast_analysis(args.query_sequence, args.target_sequences, args.output_path, args.chunk_size, args.score_type, args.score_threshold)



