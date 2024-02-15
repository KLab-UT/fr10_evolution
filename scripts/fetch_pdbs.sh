#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pdb_list_file> <target_directory>"
    exit 1
fi

pdb_list_file=$1
target_directory=$2

# Create the target directory if it doesn't exist
mkdir -p "$target_directory"

# Loop through each PDB ID in the list and download the corresponding PDB file
while IFS= read -r pdb_id; do
    # Download the PDB file
    pdb_filename="$target_directory/$pdb_id.pdb"
    echo "Downloading $pdb_id to $pdb_filename..."
    wget -O "$pdb_filename" "https://files.rcsb.org/download/$pdb_id.pdb"

    # Extract all title lines from the downloaded PDB file using awk
    pdb_titles=$(awk '/^TITLE/ {sub(/^TITLE\s*/, ""); gsub(/"/, ""); print}' "$pdb_filename")

    # Combine multiline titles into a single title, separated by space
    pdb_title_combined=$(echo "$pdb_titles" | tr -d '\n' | sed 's/TITLE[[:space:]]*2/TITLE/g')

    # Remove trailing underscores and white spaces from the combined title
    pdb_title_cleaned=$(echo "$pdb_title_combined" | sed -E 's/[[:space:]_]*$//;s/^[[:space:]_]*//')

    # If the title is empty, set it to "Untitled"
    if [ -z "$pdb_title_cleaned" ]; then
        pdb_title_cleaned="Untitled"
    fi

    # Replace spaces with underscores in the cleaned title
    pdb_title_underscored=$(echo "$pdb_title_cleaned" | tr ' ' '_')

    # Remove parentheses and their content from the cleaned title
    pdb_title_final=$(echo "$pdb_title_underscored" | sed 's/([^)]*)//g')

    # Remove double underscores, trailing underscores, and commas from the title
    pdb_title_final_cleaned=$(echo "$pdb_title_final" | sed -E 's/__+/_/g;s/[_,]+$//;s/,//g')

    # Construct the new filename with ID and cleaned title
    new_filename="$target_directory/${pdb_id}_${pdb_title_final_cleaned}.pdb"

    # Rename the downloaded PDB file
    mv "$pdb_filename" "$new_filename"
done < "$pdb_list_file"

echo "Download and renaming complete."

