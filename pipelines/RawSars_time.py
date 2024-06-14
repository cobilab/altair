import subprocess
from datetime import datetime

# Function to count sequences in a FASTA file
def count_sequences(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file if line.startswith(">"))

# Function to run the AltaiR raw command
def run_raw(input_fasta, output_fasta, min_length, max_length):
    # Construct the command as a single string
    command = f"./AltaiR raw -v -min {min_length} -max {max_length} {input_fasta} {output_fasta}"

    # Start the timer
    start_time = datetime.now()

    # Execute the command and redirect output to a file
    with open("rawoutput.fa", "w") as outfile:
        subprocess.run(command, shell=True, stdout=outfile, check=True)

    # End the timer
    end_time = datetime.now()

    # Calculate elapsed time in seconds
    elapsed_time = (end_time - start_time).total_seconds()

    return elapsed_time

# Main function to handle workflow
def main(input_fasta, output_fasta, min_length, max_length):
    total_sequences = count_sequences(input_fasta)
    elapsed_time = run_raw(input_fasta, output_fasta, min_length, max_length)
    average_time_per_sequence = elapsed_time / total_sequences if total_sequences > 0 else 0
    
    print(f"Total time taken for calculation: {elapsed_time} seconds")
    print(f"Average time per sequence: {average_time_per_sequence} seconds")

# Example usage
if __name__ == "__main__":
    main("human.fa", "sorted_output.fa", 12, 14) 

