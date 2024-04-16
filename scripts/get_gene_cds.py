import argparse
from Bio import SeqIO, Entrez

def fetch_gene_cds(gene_symbol, species_taxid):
    """
    return the accession id for the chromosome and 
    """

    search_handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene] AND {species_taxid}[Taxonomy ID]")
    gene_id_list = Entrez.read(search_handle)["IdList"]

    if not gene_id_list:
        print(f"No record found for {gene_symbol} in the Gene database for human (Taxonomy ID 9606).")
        return None
    gene_id = gene_id_list[0]

    try:
        gene_handle = Entrez.efetch(db="gene", id=gene_id, retmode="text")
        gene_record = gene_handle.read()
        print(f"gene record: {gene_record}")

    except Exception as e:
        print(f"Error fetching record for {gene_id}: {e}")
        return None
    # for feature in gene_record.features:
    #     print(f"gene feature: {feature}\nfeature type: {feature.type}\nfeature qualifiers: {feature.qualifiers}")

def main():
    parser = argparse.ArgumentParser(description="Download fasta sequence for a gene in a specified organism.")
    parser.add_argument("gene_symbol", type=str, help="Gene symbol to search for.")
    parser.add_argument("species_taxid", type=str, help="Taxonomy ID of the organism.")
    parser.add_argument("email", type=str, help="Email for NCBI Entrez.")

    args = parser.parse_args()
    Entrez.email = args.email
    fetch_gene_cds(args.gene_symbol, args.species_taxid)

if __name__ == "__main__":
    main()
