import subprocess
from datetime import datetime

# Function to count sequences in a FASTA file
def count_sequences(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file if line.startswith(">"))

# Function to run the AltaiR NCD calculation
def run_ncd(ref_sequence, fasta_file):
    # Prepare the command
    command = [
        "./AltaiR", "ncd", "-v", "--dna", "--threads", "1", 
        "-m", "12:50:0:1:0:0.9/2:10:0.9", 
        "-r", ref_sequence, fasta_file
    ]
    
    # Start the timer
    start_time = datetime.now()
    
    # Execute the command
    with open("sim-data.csv", "w") as output_file:
        subprocess.run(command, stdout=output_file, check=True)
    
    # End the timer
    end_time = datetime.now()
    
    # Calculate elapsed time in seconds
    elapsed_time = (end_time - start_time).total_seconds()
    
    return elapsed_time

# Main function to handle workflow
def main(ref_sequence, fasta_file):
    total_sequences = count_sequences(fasta_file)
    elapsed_time = run_ncd(ref_sequence, fasta_file)
    average_time_per_sequence = elapsed_time / total_sequences if total_sequences > 0 else 0
    
    print(f"Total time taken for NCD calculation: {elapsed_time} seconds")
    print(f"Average time per sequence: {average_time_per_sequence} seconds")

# Example usage
if __name__ == "__main__":
    main("ORIGINAL.fa", "ALL.fa")
