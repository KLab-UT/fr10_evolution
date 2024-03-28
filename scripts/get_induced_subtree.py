import json
import sys
import argparse
import subprocess
import re

def extract_ott_ids(json_file):
    """
    Extracts OTT IDs from the TNRS nameset JSON output.

    Args:
        json_file (str): Path to the TNRS nameset JSON file.

    Returns:
        list: List of OTT IDs.
    """
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
            ott_ids = [entry.get('ottId') for entry in data.get('names', []) if entry.get('ottId') is not None]
            return ott_ids
    except FileNotFoundError:
        print(f"Error: File '{json_file}' not found.")
    except json.JSONDecodeError:
        print(f"Error: Unable to decode JSON in file '{json_file}'. Please ensure it's valid JSON.")

def get_induced_subtree(id_list, induced_json_output, output_treefile):
    """
    Retrieves the induced subtree from Open Tree of Life using the provided OTT IDs.

    Args:
        id_list (list): List of OTT IDs.
        induced_json_output (str): Path to store the induced subtree JSON output.
        output_treefile (str): Path to store the fixed tree.

    Returns:
        None
    """
    try:
        # Convert id_list to a JSON array
        id_list_json = json.dumps(id_list)

        # Get the induced subtree
        command = f"curl -X POST https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree -H 'content-type: application/json' -d '{{\"ott_ids\":{id_list_json}}}' > {induced_json_output}"
        subprocess.run(command, shell=True, check=True)

        # Load induced subtree JSON and extract the newick tree if it exists
        with open(induced_json_output, 'r') as json_file:
            try:
                tree_data = json.load(json_file)
                tree = tree_data.get("newick")

                # Check if 'newick' key exists before manipulation
                if tree is not None:
                    # Remove anything within single quotes from the tree
                    fixed_tree = re.sub(r"'[^']*'", '', tree)

                    # Remove '_ott' from the tree and write to the output file
                    fixed_tree = re.sub(r'_ott\d+', '', fixed_tree)
                    with open(output_treefile, "w") as output_tree:
                        output_tree.write(fixed_tree)
                else:
                    print("Error: 'newick' key not found in the induced subtree JSON.")
                    # Print the content of the induced JSON for inspection
                    print("Induced JSON content:", tree_data)

            except json.JSONDecodeError:
                print("Error: Unable to decode JSON in induced subtree JSON file.")

    except subprocess.CalledProcessError as e:
        print("Error: Unable to fetch induced subtree from Open Tree of Life.")
        # Print the error message and return code
        print("Error Message:", e.output.decode())
        print("Return Code:", e.returncode)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"Create induced subtree from OTOL tnrs nameset json output\nSee: https://tree.opentreeoflife.org/curator/tnrs/")
    parser.add_argument("input_json_file", type=str, help="Path to the TNRS nameset JSON file")
    parser.add_argument("induced_json_output", type=str, help="Path to store the induced subtree JSON output")
    parser.add_argument("output_file", type=str, help="Path to store the fixed subtree")
    
    args = parser.parse_args()

    ott_ids = extract_ott_ids(args.input_json_file)
    get_induced_subtree(ott_ids, args.induced_json_output, args.output_file)
