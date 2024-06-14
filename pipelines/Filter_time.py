import subprocess
from datetime import datetime

# Calculate the number of sequences
def count_sequences(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file if line.startswith(">"))

# Call the external AltaiR filter command
def filter_sequences(min_val, max_val):
    start_time = datetime.now()
    command = f"cat SARS-CoV-2.fa | ./AltaiR filter -a ACGT -min {min_val} -max {max_val} -i '|2019|' -i '|2020|' -i '|2021|' -i '|2022|' " \
              f"| ./AltaiR filter -p '|2019' -p '|2020' -p '|2021' -p '|2022' -p '2019|' -p '2020|' -p '2021|' -p '2022|' > FIL-SARS-CoV-2.fa"
    subprocess.run(command, shell=True, check=True)
    end_time = datetime.now()
    elapsed_time = (end_time - start_time).total_seconds()

    return elapsed_time

# Main function
def main(min_val, max_val):
    total_sequences = count_sequences("SARS-CoV-2.fa")
    elapsed_time = filter_sequences(min_val, max_val)
    average_time_per_sequence = elapsed_time / total_sequences if total_sequences else 0

    print(f"Total time taken for filtering: {elapsed_time} seconds")
    print(f"Average time per sequence: {average_time_per_sequence} seconds")

# Example usage
main(29885, 29921)

#Total time taken for filtering: 598.435644 seconds
#Average time per sequence: 0.00038907586592505666 seconds
#Total time taken for NC calculation: 128.25567 seconds
#Average time per sequence: 0.005011161600375088 seconds
#Total time taken for frequency calculation: 3.701271 seconds
#Average time per sequence: 0.0001446147925294991 seconds
#Total time taken for NCD calculation: 2294.131466 seconds
#Average time per sequence: 0.22941314659999998 seconds