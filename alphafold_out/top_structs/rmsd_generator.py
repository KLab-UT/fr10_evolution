from pymol import cmd

cmd.load('fr10_r_syl.pdb')
cmd.load('apoA-II_x_trop.pdb')
cmd.show('spheres','all')
cmd.png('output_image.png', width=800, height=800, dpi=200)
cmd.align('fr10_r_syl', 'apoA-II_x_trop')
