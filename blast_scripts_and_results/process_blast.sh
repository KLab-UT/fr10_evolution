#Takes a textfile of line seperated blast outputs in outfmt 10 w/ seq as last value for each hit, extracts hit seqs into temp.fasta file then splits each into indv. fastas w/ hit seq and query seq 

infile=$1
query=$2
min_length=$3


python3 get_seqs.py $infile $infile.temp.fasta $min_length
bash combine.sh $infile.temp.fasta $query ${infile}_tempdir

#ex usage : bash process_blast.sh tblastx_exons_results.txt fr10.faa 10
