
#download the genomes
module load blast/2.11.0
#cd ~/fr10_evolution
#bash /uufs/chpc.utah.edu/common/home/u6052680/fr10_evolution/get_refs.sh -d /scratch/general/nfs1/utu_4310/fr10_evolution_wd/reference_genomes -l fixed_ref_genbank.txt 
#gunzip -k /scratch/general/nfs1/utu_4310/fr10_evolution_wd/reference_genomes/*.gz
#mkdir /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs
#mv /scratch/general/nfs1/utu_4310/fr10_evolution_wd/reference_genomes/*.fna /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs

#make blast db for each genome

#for genome in /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs/*.fna; do makeblastdb -in $genome -dbtype nucl -out $genome.db
#done

makeblast db -in nanpar.fasta -dbtype nucl -out nanpar.fasta.db
blastn -db nanpar.fasta.db -query fr10.fasta -evalue 0.001 -outfmt fmt "10 stitle qseqid sseqid evalue sstart send sseq" -out nanpar_blastn_results.txt

#blast each genome for fr10
#for genome in /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs/*.fna; do blastn -query  /uufs/chpc.utah.edu/common/home/u6052680/fr10_evolution/fr10_exin.fasta -db "$genome.db" -evalue 0.001 -outfmt "10 stitle qseqid sseqid evalue sstart send sseq" -out "$genome.exin.results.txt"
#done
#mv /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs/*.db.results.txt /scratch/general/nfs1/utu_4310/fr10_evolution_wd/blast_outputs
