##Need to download fr10/scratch/general/nfs1/utu_4310/fr10_evolution_wd/reference_genome.fasta
##Need to make a directory for all the dbs

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

#blast each genome for fr10
for genome in /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs/*.fna; do blastn -query  /uufs/chpc.utah.edu/common/home/u6052680/fr10_evolution/fr10.fasta -db "$genome.db" -evalue 0.001 -outfmt "5" -out "$genome.db.results.xml"
done
mv /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbs/*.db.results.txt /scratch/general/nfs1/utu_4310/fr10_evolution_wd/blast_outputs
