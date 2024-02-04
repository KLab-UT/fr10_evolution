import re
import argparse

def check_parentheses_balance(treefile):
    left_parentheses_count = treefile.count('(')
    right_parentheses_count = treefile.count(')')

    return left_parentheses_count, right_parentheses_count

def replace_ids_with_descriptions(treefile):
    # Split the treefile into lines
    lines = treefile.split(',')

    # Define a dictionary to store sequence IDs and their descriptions
    id_description_map = {}

    # Extract sequence IDs and descriptions from each line
    for line in lines:
        match = re.match(r'(\S+)\s+(\S.+)', line)
        if match:
            sequence_id = match.group(1)
            description = match.group(2)
            id_description_map[sequence_id] = description

    # Replace sequence IDs with descriptions
    for sequence_id, description in id_description_map.items():
        treefile = treefile.replace(sequence_id, description)

    return treefile

def main():
    parser = argparse.ArgumentParser(description='Replace sequence IDs with descriptions in a treefile.')
    parser.add_argument('input_file', type=str, help='Path to the input treefile')
    parser.add_argument('output_file', type=str, help='Path to the output treefile')

    args = parser.parse_args()

    try:
        with open(args.input_file, 'r') as f:
            treefile_content = f.read()

        left_count, right_count = check_parentheses_balance(treefile_content)

        if left_count != right_count:
            print("Error: The number of left and right parentheses in the Newick string is not equal.")
            return

        updated_treefile = replace_ids_with_descriptions(treefile_content)

        with open(args.output_file, 'w') as f:
            f.write(updated_treefile)

        print(f"Replacement successful. Updated treefile saved to {args.output_file}")

    except FileNotFoundError:
        print("File not found. Please provide a valid file path.")

if __name__ == "__main__":
    main()

