#!/bin/bash

measure_command() {
    local name=$1
    local command=$2
    echo "Running $name test..."
    
    # Start time
    start=$(date +%s.%N)
    
    # Run command and capture output
    output=$(eval "$command" 2>&1)
    exit_code=$?
    
    # End time
    end=$(date +%s.%N)
    
    # Calculate runtime
    runtime=$(echo "$end - $start" | bc)
    
    # Get peak memory usage (in KB, convert to MB)
    peak_memory=$(echo "$output" | grep "Maximum resident set size" | awk '{print $6/1024}')
    
    echo "$name - Total time: $runtime s, Peak RAM: $peak_memory MB"
    
    # Save results to JSON file
    echo "{\"total_time\": $runtime, \"peak_ram\": $peak_memory}" > "${name}_result.json"
    
    if [ $exit_code -ne 0 ]; then
        echo "Error running $name test. Output:"
        echo "$output"
    fi
}

# Filtering
measure_command "Filtering" "/usr/bin/time -v bash -c 'cat SARS-CoV-2.fa | ./AltaiR filter -a ACGT -min 29885 -max 29921 -i \"|2019|\" -i \"|2020|\" -i \"|2021|\" -i \"|2022|\" | ./AltaiR filter -p \"|2019\" -p \"|2020\" -p \"|2021\" -p \"|2022\" -p \"2019|\" -p \"2020|\" -p \"2021|\" -p \"2022|\" > FIL-SARS-CoV-2.fa'"

# NC
measure_command "NC" "/usr/bin/time -v ./AltaiR nc -v --dna --threads 1 -m 2:1:0:1:0:0.8/0:0:0 -m 3:1:0:0:0:0.9/0:0:0 FIL-SARS-CoV-2.fa > comp-data.csv"

# Frequency
measure_command "Frequency" "/usr/bin/time -v ./AltaiR frequency -v -a ACTG < sorted_output.fa > freq-data.csv"

# NCD
measure_command "NCD" "/usr/bin/time -v ./AltaiR ncd -v --dna --threads 1 -m 12:50:0:1:0:0.9/2:10:0.9 -r ORIGINAL.fa ALL.fa > sim-data.csv"

# RAW
measure_command "RAW" "/usr/bin/time -v ./AltaiR raw -v -min 12 -max 14 humanGenomeAndTranscriptome.fa sorted_output.fa > rawoutput.fa"