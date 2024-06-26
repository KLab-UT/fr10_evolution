#Reload trees 
git checkout HEAD -- ../apo-fr10_alignments/

fr10_tree_files=(
  "../alignments/apo-fr10_alignments/ApoA-II_aligned.plottree"
  "../alignments/apo-fr10_alignments/ApoA-V_aligned.plottree"
  "../alignments/apo-fr10_alignments/ApoC-IV_aligned.plottree"
  "../alignments/apo-fr10_alignments/ApoA-IV_aligned.plottree"
  "../alignments/apo-fr10_alignments/ApoC-III_aligned.plottree"
  "../alignments/apo-fr10_alignments/ApoC-I_aligned.plottree"
  "../alignments/apo-fr10_alignments/ApoA-I_aligned.plottree"
  "../alignments/apo-fr10_alignments/ApoC-II_aligned.plottree"
  "../alignments/apo-fr10_alignments/АроЕ_aligned.plottree"
)

drp10_tree_files=(
"../alignments/apo-drp10_alignments/ApoA-II_aligned.plottree"
"../alignments/apo-drp10_alignments/ApoA-V_aligned.plottree"
"../alignments/apo-drp10_alignments/ApoC-IV_aligned.plottree"
"../alignments/apo-drp10_alignments/ApoA-IV_aligned.plottree"
"../alignments/apo-drp10_alignments/ApoC-III_aligned.plottree"
"../alignments/apo-drp10_alignments/ApoC-I_aligned.plottree"
"../alignments/apo-drp10_alignments/ApoA-I_aligned.plottree"
"../alignments/apo-drp10_alignments/ApoC-II_aligned.plottree"
"../alignments/apo-drp10_alignments/АроЕ_aligned.plottree"
)


# Setup output files
header="Gene,Homo_sapiens,Pan_troglodytes,gorilla_gorilla,Pongo_abelii,Macaca_mulatta,Papio_anubis,Mus_musculus,Rattus_norvegicus,Oryctolagus_cuniculus,Bos_taurus,Ovis_aries,Capra_hircus,Sus_scrofa,Equus_caballus,Felis_catus,lupus_familiaris,Monodelphis_domestica,Ornithorhynchus_anatini,Gallus_gallus,Meleagris_gallopavo,Anolis_carolinensis,Chrysemys_picta,Takifugu_rubripes,Oryzias_latipes,Danio_rerio,Lepisosteus_oculatus,Latimeria_chalumnae,Rhincodon_typus,Callorhinchus_milii"
echo $header > normalized_distances_fr10.csv
echo $header > normalized_distances_drp10.csv


# Call R script to calculate and append normalized distances for each treefile and append the distances for that tree to the normalized_distances output"
for tree_file in "${fr10_tree_files[@]}"; do
  echo "$tree_file"
  echo $header > normalized_distances_fr10.csv
  Rscript plot_phylo_dist.r "$tree_file" normalized_distances_fr10.csv Lithobates_sylvaticus;
done


















