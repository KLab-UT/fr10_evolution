# Importing molecular structure from pymol
from pymol import cmd
# Comparing pbd files to get RMSD values
def compare_rmsd(selection1, selection2):
    """
    Calculate RMSD between two selections.
    """
    # Align the two selections and return the RMSD value
    rmsd = cmd.align(selection1, selection2)[0]
    return rmsd

# List of PDB files
pdb_files = [
    "fr10_r_den.pdb",
    "apoA-I_x_trop.pdb",
    "apoA-II_x_trop.pdb",
    "apoA-IV_x_trop.pdb",
    "apoA-V_x_trop.pdb",
    "apoC-I_x_trop.pdb",
    "apoC-II_x_trop.pdb",
    "apoE_x_trop.pdb",
    "drp10_x_laev.pdb",
    "fr10_r_syl.pdb",
]

outfile = open('rmsd_values.txt', 'w')

# Function to compare each file to each other
def compare_all_rmsd():
    for i, file1 in enumerate(pdb_files):
        for j, file2 in enumerate(pdb_files):
            if i < j:
                cmd.reinitialize()
                cmd.load(file1, "prot1")
                cmd.load(file2, "prot2")
                rmsd = compare_rmsd("prot1", "prot2")
                outfile.write("RMSD between '{}' and '{}': {:.3f} Å".format(file1, file2, rmsd))

# Execute comparison
compare_all_rmsd()

outfile.close()
