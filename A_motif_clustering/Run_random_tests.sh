#!/bin/bash -l

module load bioinfo-tools

# Paths
path_hammock="/home/juliaa/projekt/Hammock/dist/Hammock.jar" # Path to hammock script
path_to_random_files="/home/juliaa/projekt/Hammock/input/random_test/random_inputs" # Path to randomised fasta files
path_to_sequences_input="/home/juliaa/projekt/Hammock/input/random_test/output_no_duplicates.fa" # Path to original fasta files without random duplicates
output_duplicates_folder="/home/juliaa/projekt/Hammock/results_new/res_random/duplicates" # Path to original fasta files that had random duplicates
output_no_duplicates_folder="/home/juliaa/projekt/Hammock/input/random_test/input_no_duplicates" #
path_final="/home/juliaa/projekt/Hammock/results_new/res_random/final_outputs" # Path to final outputs
path_initial_clusters="/home/juliaa/projekt/Hammock/results_new/res_random/initial_clusters" # Path to inital cluster cores created by the randomised files 



# Process files in the random input directory
for file in "$path_to_random_files"/*; do
  filename=$(basename "$file" .fa)  # Get the base name of the file without extension

  # Generate the output file paths
  duplicates_file="$output_duplicates_folder/duplicates_${filename}.fa"
  no_duplicates_file="$output_no_duplicates_folder/no_duplicates_${filename}.fa"

  # Remove duplicate sequences 
  grep -v "^>" "$file" | sort | uniq > sequences_file.txt

  # Compare $path_to_sequences_input against $file
  awk 'BEGIN {while (getline < "sequences_file.txt") seen[$0] = 1}
     /^>/ {h=$0; getline s;
           if (s in seen) {print h > "'"$duplicates_file"'"; print s > "'"$duplicates_file"'"}
           else {print h > "'"$no_duplicates_file"'"; print s > "'"$no_duplicates_file"'"}}' "$path_to_sequences_input"
  rm sequences_input.txt
	
  java -jar "$path_hammock"  full -i "$file" —-gap_penalty 0 —-order random —seed 42 --max_shift 2  -d "$path_initial_clusters/i_${filename}"
  java -jar "$path_hammock" cluster -i "$path_initial_clusters/i_${filename}/final_clusters.tsv" -as "$no_duplicates_file" —-gap_penalty 0 —-order random —seed 42 --max_shift 2  -d "$path_final/f_${filename}"

done
