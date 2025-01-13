#!/bin/bash

# Script to get the random clusters created from the ... script

# Define source and target directories
source_dir="/home/juliaa/projekt/Hammock/results_new/res_random/final_outputs" 
target_dir="/home/juliaa/projekt/Hammock/results_new/res_random/final_clusters"

    # Define the source file path
    source_file="$dir/final_clusters.tsv"

    # Define the target file path
    target_file="$target_dir/clusters_${directory_name}.tsv"

    # Move and rename the file if it exists
    if [ -f "$source_file" ]; then
      mv "$source_file" "$target_file"
    fi
  fi
done

