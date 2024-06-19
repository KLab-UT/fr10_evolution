from pymol import cmd

def compare_rmsd(selection1, selection2):
    """
    Calculate RMSD between two selections.
    """
    rmsd = cmd.align(selection1, selection2)[0]
    return rmsd

# List of PDB files
pdb_files = [
    "fr10_r_den.pdb",
    "apoA-I_x_trop (1).pdb",
    "apoA-II_x_trop (1).pdb",
    "apoA-IV_x_trop (1).pdb",
    "apoA-V_x_trop (1).pdb",
    "apoC-I_x_trop (1).pdb",
    "apoC-II_x_trop (1).pdb",
    "apoE_x_trop (1).pdb",
    "drp10_x_laev (1).pdb",
    "fr10_a_wuy.pdb",
    "fr10_l_clam.pdb",
    "fr10_l_pip.pdb",
    "fr10_o_tor.pdb",
    "fr10_p_ads.pdb",
    "fr10_p_mega.pdb",
    "fr10_r_cat.pdb",
    "fr10_r_kuk.pdb",
    "fr10_r_mus.pdb",
    "fr10_r_syl.pdb",
    "fr10_r_temp.pdb",
    "ls12_m_oct.pdb"
]

# Function to compare each file to each other
def compare_all_rmsd():
    for i, file1 in enumerate(pdb_files):
        for j, file2 in enumerate(pdb_files):
            if i < j:
                cmd.reinitialize()
                cmd.load(file1, "prot1")
                cmd.load(file2, "prot2")
                rmsd = compare_rmsd("prot1", "prot2")
                print("RMSD between '{}' and '{}': {:.3f} Å".format(file1, file2, rmsd))

# Execute comparison
compare_all_rmsd()

compare_rmsd(selection1, selection2)



#New Code 