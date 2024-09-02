import subprocess
from datetime import datetime
import os
import sys

def count_sequences(filename):
    print(f"Counting sequences in {filename}...")
    try:
        with open(filename, 'r') as file:
            count = sum(1 for line in file if line.startswith(">"))
        print(f"Found {count} sequences.")
        return count
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        sys.exit(1)

def run_command(command):
    print(f"Running command: {command}")
    try:
        subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)

def filter_sequences_altair(min_val, max_val):
    print("Starting AltaiR filtering...")
    start_time = datetime.now()
    altair_path = "./AltaiR"
    if not os.path.exists(altair_path):
        print(f"Error: AltaiR executable not found at {altair_path}")
        sys.exit(1)
    command = f"cat SARS-CoV-2.fa | {altair_path} filter -a ACGT -min {min_val} -max {max_val} -i '|2019|' -i '|2020|' -i '|2021|' -i '|2022|' " \
              f"| {altair_path} filter -p '|2019' -p '|2020' -p '|2021' -p '|2022' -p '2019|' -p '2020|' -p '2021|' -p '2022|' > FIL-SARS-CoV-2.fa"
    run_command(command)
    end_time = datetime.now()
    return (end_time - start_time).total_seconds()

def filter_sequences_grep(min_val, max_val):
    print("Starting grep filtering...")
    start_time = datetime.now()
    command = f"""
    awk -v min={min_val} -v max={max_val} '
    BEGIN {{RS=">";FS="\\n"}} 
    NR>1 {{
        sequence=$2;
        for(i=3;i<=NF;i++) sequence=sequence $i;
        if(length(sequence) >= min && length(sequence) <= max && sequence ~ /^[ACGT]+$/ && 
           ($1 ~ /2019/ || $1 ~ /2020/ || $1 ~ /2021/ || $1 ~ /2022/) && 
           $1 !~ /\|2019\|/ && $1 !~ /\|2020\|/ && $1 !~ /\|2021\|/ && $1 !~ /\|2022\|/)
            print ">"$0
    }}' SARS-CoV-2.fa > GREP-FIL-SARS-CoV-2.fa
    """
    run_command(command)
    end_time = datetime.now()
    return (end_time - start_time).total_seconds()

def main(min_val, max_val):
    total_sequences = count_sequences("SARS-CoV-2.fa")

    # AltaiR filtering
    altair_time = filter_sequences_altair(min_val, max_val)
    altair_avg_time = altair_time / total_sequences if total_sequences else 0

    print(f"\nAltaiR filtering:")
    print(f"Total time: {altair_time:.6f} seconds")
    print(f"Average time per sequence: {altair_avg_time:.6f} seconds")

    # Grep/Awk filtering
    grep_time = filter_sequences_grep(min_val, max_val)
    grep_avg_time = grep_time / total_sequences if total_sequences else 0

    print(f"\nGrep/Awk filtering:")
    print(f"Total time: {grep_time:.6f} seconds")
    print(f"Average time per sequence: {grep_avg_time:.6f} seconds")

    # Comparison
    speedup = grep_time / altair_time if altair_time > 0 else float('inf')
    print(f"\nAltaiR is {speedup:.2f}x faster than Grep/Awk for this filtering task.")

if __name__ == "__main__":
    main(29885, 29921)


# Grep/Awk filtering:
# Total time: 720.381350 seconds
# Average time per sequence: 0.000816 seconds

# AltaiR is 5.55x faster than Grep/Awk for this filtering task.