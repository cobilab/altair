import sys
import random

def process_sequence(seq):
    """
    Process a sequence: convert to uppercase and replace any character not in 'ACGT' with a random choice from 'ACGT'.
    """
    return ''.join(random.choice('ACGT') if c not in 'ACGT' else c.upper() for c in seq)

def count_bases(sequences):
    """
    Count the number of each base in the sequences.
    """
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for seq in sequences:
        for base in seq:
            if base in counts:
                counts[base] += 1
    return counts

def write_processed_fasta(sequences, headers, output_file_path):
    """
    Write processed sequences to a new fasta file.
    """
    with open(output_file_path, 'w') as file:
        for header, sequence in zip(headers, sequences):
            file.write(header + '\n')
            file.write(sequence + '\n')

def process_fasta(input_file_path, output_file_path):
    """
    Process a multi-fasta file and write to a new file.
    """
    sequences = []
    headers = []
    current_sequence_parts = []

    with open(input_file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if current_sequence_parts:
                    sequences.append(process_sequence(''.join(current_sequence_parts)))
                    current_sequence_parts = []
                headers.append(line.strip())
                print("Processing:", line.strip())  # Progress tracking
            else:
                current_sequence_parts.append(line.strip())

        if current_sequence_parts:
            sequences.append(process_sequence(''.join(current_sequence_parts)))

    write_processed_fasta(sequences, headers, output_file_path)

    base_counts = count_bases(sequences)
    return base_counts, ''.join(sorted(base_counts.keys()))

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file_path> <output_file_path>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    base_counts, alphabet = process_fasta(input_file_path, output_file_path)
    print("Base Counts:", base_counts)
    print("Alphabet:", alphabet)

if __name__ == "__main__":
    main()
