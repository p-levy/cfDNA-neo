#!/bin/bash

set -uex

# Function to display usage
usage() {
  echo "Usage: $0 -i INPUT_FILE -r R_SCRIPT"
  exit 1
}

# Function to run R script for a single line
run_rscript() {
  local line="$1"
  # Remove carriage return characters from the line
  line=$(echo "$line" | tr -d '\r')
  # Read the fields from the CSV line
  IFS=',' read -r field1 field2 field3 field4 field5 field6 <<< "$line"
  # Run the R script with the extracted fields
  Rscript "$R_SCRIPT" "$field1" "$field2" "$field3" "$field4" "$field5" "$field6"
}

# Processing function
run() {
  tail -n +2 "$INPUT_FILE" | while IFS= read -r line || [[ -n "$line" ]]; do
    run_rscript "$line"
  done
}

# Parse arguments
while getopts "i:r:" opt; do
  case $opt in
    i)
      INPUT_FILE="$OPTARG"
      ;;
    r)
      R_SCRIPT="$OPTARG"
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

# Run
run