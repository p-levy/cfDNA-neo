#!/bin/bash

set -e

# Function to display usage
usage() {
  echo "Usage: $0 -i INPUT_FILE -r R_SCRIPT [-p NUM_CPUS]"
  echo "INPUT_FILE is a csv with the field: patient,variants_counts_path,outdir,segs_cfDNA,segs_FrTu,bed_exome,pipeline_excel,immunogenic,purity_FrTu,purity_cfDNA"
  exit 1
}



# Function to run R script for a single line
run_rscript() {
  local line="$1"
  # Remove carriage return characters from the line
  line=$(echo "$line" | tr -d '\r')
  # Read the fields from the CSV line
  IFS=',' read -r field1 field2 field3 field4 field5 field6 field7 field8 field9 field10 <<< "$line"
  # Run the R script with the extracted fields
  Rscript "$R_SCRIPT" "$field1" "$field2" "$field3" "$field4" "$field5" "$field6" "$field7" "$field8" "$field9" "$field10"
}

# Sequential processing function
run_sequential() {
  tail -n +2 "$INPUT_FILE" | while IFS= read -r line || [[ -n "$line" ]]; do
    run_rscript "$line"
  done
}

# Parallel processing function
run_parallel() {
  local num_cpus="$1"
  export -f run_rscript
  export R_SCRIPT
  tail -n +2 "$INPUT_FILE" | parallel -j "$num_cpus" run_rscript
}

# Parse arguments
while getopts "i:r:p:" opt; do
  case $opt in
    i)
      INPUT_FILE="$OPTARG"
      ;;
    r)
      R_SCRIPT="$OPTARG"
      ;;
    p)
      MODE="parallel"
      NUM_CPUS="$OPTARG"
      ;;
    *)
      usage
      ;;
  esac
done

# Check for mandatory arguments
if [[ -z "$INPUT_FILE" || -z "$R_SCRIPT" ]]; then
  usage
fi

# Default to sequential mode if -p is not provided
MODE="${MODE:-sequential}"
NUM_CPUS="${NUM_CPUS:-1}"

# Check if input file exists
if [[ ! -f $INPUT_FILE ]]; then
  echo "Error: Input file '$INPUT_FILE' not found."
  exit 1
fi

# Check if R script exists
if [[ ! -f $R_SCRIPT ]]; then
  echo "Error: R script '$R_SCRIPT' not found."
  exit 1
fi

# Run the appropriate mode
if [[ $MODE == "sequential" ]]; then
  run_sequential
else
  run_parallel "$NUM_CPUS"
fi