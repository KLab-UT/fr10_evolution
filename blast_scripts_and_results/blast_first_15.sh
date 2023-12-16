'''script for searching for primer sequences w/ blastn -task "blastn-short", loops through each genome given and search for provided primer sequences'''

module load blast
genomes=/uufs/chpc.utah.edu/common/home/u6052680/fr10_evolution/hit
primers=/uufs/chpc.utah.edu/common/home/u6052680/fr10_evolution/primers.fna

while IFS= read -r line; do
  # Perform actions on each line here
  blastn -db  "/scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs/$line.db" -query "/uufs/chpc.utah.edu/common/home/u6052680/fr10_evolution/fr10_first_15.fasta" -evalue 0.001 -outfmt "10 stitle qtitle qseqid sseqid evalue sstart send sseq length pident mismatch frames" -out "/scratch/general/nfs1/utu_4310/fr10_evolution_wd/blast_outputs/$line.first.results.txt"
done < "$genomes"
