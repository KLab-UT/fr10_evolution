alignments_dir=$1
query_species=$2

#get alignments to make plots
for alignment in "$alignments_dir"/*aligned.fa; do
    echo "Making plots for: $alignment"

    # Create variables for .plottree and _distmatrix.csv for formatted tree file(species only) and distance matrix
    plottree_file="${alignment%.fa}.plottree"
    distmatrix_file="${alignment%.fa}_distmatrix.csv"
    echo "Plottree file: $plottree_file"
    echo "Distmatrix file: $distmatrix_file"

    # Format tree descriptions to only contain the species name for each sequence
    python3 scripts/get_tree_species.py "${alignment}.contree" "$alignment" "$plottree_file"

    # Call R script to get distance matrix
    #Rscript plot_phylo_dist.r $plottree_file $distmatrix_file
done


##### USAGE FOR REGENERATING apo-fr10 and apo-drp10 plot trees
# bash scripts/make_plots.sh alignments/apo-drp10_alignments Xenopus_laevis
# bash scripts/make_plots.sh alignments/apo-fr10_alignments Lithobates_sylvaticus
