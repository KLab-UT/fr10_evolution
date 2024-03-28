import json
import argparse
from Bio import Phylo

def switch_tree_names(species_map, treefile, outfile):
    """
    Switch the names in a phylogenetic tree file based on the provided species map.

    Parameters:
    - species_map (str): Path to the JSON file containing the species mapping.
    - treefile (str): Path to the input tree file.
    - outfile (str): Path to the output file to save the modified tree.

    Returns:
    None
    """
    # Load species map from JSON file
    with open(species_map, "r") as file:
        mapping = json.load(file)

    # Read the input tree file
    with open(treefile, "r") as tree_file:
        tree = tree_file.read()

    # Replace species names in the tree
    for species in mapping:
        if species in tree:
            print(f"Replacing {species} with {mapping[species]}")
            tree = tree.replace(species, mapping[species])

    # Save the modified tree to the output file
    with open(outfile, "w") as out_file:
        out_file.write(tree)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Switch names in a phylogenetic tree based on a species map.")
    parser.add_argument("species_map", help="Path to the JSON file containing the species mapping.")
    parser.add_argument("treefile", help="Path to the input tree file.")
    parser.add_argument("outfile", help="Path to the output file to save the modified tree.")
    
    # Parse command line arguments
    args = parser.parse_args()

    # Call the switch_tree_names function with provided arguments
    switch_tree_names(args.species_map, args.treefile, args.outfile)

if __name__ == "__main__":
    main()


    

    

