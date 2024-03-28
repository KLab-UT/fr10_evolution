"""parse_blast_xml.py
processes complete sequences from webblast xml results intedted to be used with blastn and tblastn for creating nucleuotide alignments
"""
import argparse
import xml.etree.ElementTree as ET
import subprocess
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class BlastHit:
    """
    Manage blast hits from ET element
    """
    def __init__(self, hit_element, blast_type):
        self.accession = hit_element.find("Hit_accession").text
        self.definition = hit_element.find("Hit_def").text
        self.length = int(hit_element.find("Hit_len").text)

        if "TSA" in self.definition:
            self.species = " ".join(self.definition.split()[:3])
        else:
            self.species = " ".join(self.definition.split()[:2])

        print(f"species for {self.accession}: {self.species}")
        self.hsps = []

        # Iterate over Hsp elements and create Hsp objects
        for hsp_element in hit_element.findall("Hit_hsps/Hsp"):
            if blast_type == "blastn":
                hsp = self.create_blastn_hsp_object(hsp_element)
            elif blast_type == "tblastn" or blast_type == "tblastx":
                hsp = self.create_tblastn_hsp_object(hsp_element)
            else:
                raise ValueError(f"Unsupported blast type: {blast_type}")
            self.hsps.append(hsp)

    def create_blastn_hsp_object(self, hsp_element):
        """
        Create Hsp object for blastn hits from Hsp XML element
        """
        hit_start=int(hsp_element.find("Hsp_hit-from").text)
        hit_end=int(hsp_element.find("Hsp_hit-to").text)
        query_start=int(hsp_element.find("Hsp_query-from").text)
        query_end=(hsp_element.find("Hsp_query-to").text)
        hseq = fetch_sequence(self.accession, hit_start, hit_end)
        return Hsp(
            bit_score=float(hsp_element.find("Hsp_bit-score").text),
            score=int(hsp_element.find("Hsp_score").text),
            evalue=float(hsp_element.find("Hsp_evalue").text),
            query_start=query_start,
            query_end=query_end,
            hit_start=hit_start,
            hit_end=hit_end,
            query_frame=int(hsp_element.find("Hsp_query-frame").text),
            hit_frame=int(hsp_element.find("Hsp_hit-frame").text),
            identity=int(hsp_element.find("Hsp_identity").text),
            positive=int(hsp_element.find("Hsp_positive").text),
            gaps=int(hsp_element.find("Hsp_gaps").text),
            align_len=int(hsp_element.find("Hsp_align-len").text),
            qseq=hsp_element.find("Hsp_qseq").text,
            hseq=hseq,
            midline=hsp_element.find("Hsp_midline").text,
        )
   
    def create_tblastn_hsp_object(self, hsp_element):
        hit_start=int(hsp_element.find("Hsp_hit-from").text)
        hit_end=int(hsp_element.find("Hsp_hit-to").text)
        query_start=3*int(hsp_element.find("Hsp_query-from").text)
        query_end=3*int(hsp_element.find("Hsp_query-to").text)
        hseq = fetch_sequence(self.accession, hit_start, hit_end)
        return Hsp(
            bit_score=float(hsp_element.find("Hsp_bit-score").text),
            score=int(hsp_element.find("Hsp_score").text),
            evalue=float(hsp_element.find("Hsp_evalue").text),
            query_start=query_start,
            query_end=query_end,
            hit_start=hit_start,
            hit_end=hit_end,
            query_frame=int(hsp_element.find("Hsp_query-frame").text),
            hit_frame=int(hsp_element.find("Hsp_hit-frame").text),
            identity=int(hsp_element.find("Hsp_identity").text),
            positive=int(hsp_element.find("Hsp_positive").text),
            gaps=int(hsp_element.find("Hsp_gaps").text),
            align_len=int(hsp_element.find("Hsp_align-len").text),
            qseq=hsp_element.find("Hsp_qseq").text,
            hseq=hseq,
            midline=hsp_element.find("Hsp_midline").text,
        )


    def get_fasta_sequences(self, query_length):
        """
        Get FASTA-formatted sequences for each HSP's hit sequence
        """
        fasta_sequences = []
        for hsp in self.hsps:
            fasta_sequences.append(hsp.get_fasta_sequence(self.accession, self.species, query_length))
        return fasta_sequences

        
class Hsp:
    """
    Manage HSP information
    An HSP is a high scoring pair representing a specific aligngment within a BLAST hit of a sequence
    """
    def __init__(self, bit_score, score, evalue, query_start, query_end,
                 hit_start, hit_end, query_frame, hit_frame, identity,
                 positive, gaps, align_len, qseq, hseq, midline):
        self.bit_score = bit_score
        self.score = score
        self.evalue = evalue
        self.query_start = query_start
        self.query_end = query_end
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.query_frame = query_frame
        self.hit_frame = hit_frame
        self.identity = identity
        self.positive = positive
        self.gaps = gaps
        self.align_len = align_len
        self.qseq = qseq
        self.hseq = hseq
        self.midline = midline

    def get_fasta_sequence(self, hit_accession, species, query_length):
        """
        Get FASTA-formatted sequence for HSP hits positioned at where it hit in the query
        """
        seq = "-" * (self.query_start - 1) + self.hseq 
        seq += ("-" * (query_length-len(seq)))

        print(f"hit sequence: {seq} /n length: {len(seq)}")

        header = f">{hit_accession}"
        description = f"{self.hit_start}-{self.hit_end} | {species}"
        seq_record = SeqRecord(Seq(seq), id=header, description=description)
        return seq_record



def get_sequence_length(query_fasta):
    """
    Get the length of the sequence in a FASTA file.

    Parameters:
    - fasta_file (str): Path to the FASTA file.

    Returns:
    - int: Length of the sequence.
    """
    with open(query_fasta, "r") as file:
        lines = file.readlines()

    # Remove newline characters and check if the file is not empty
    lines = [line.strip() for line in lines if line.strip()]
    if not lines:
        raise ValueError("FASTA file is empty.")

    # Assume the first line is the header, and concatenate the rest as the sequence
    sequence = "".join(lines[1:])

    return len(sequence)

def fetch_sequence(accession_number, start, end):
    """
    Fetch a specified range of sequence from NCBI for a given accession number.

    Parameters:
    - email (str): Your email address for NCBI.
    - accession_number (str): Accession number to fetch.
    - start (int): Start position of the sequence.
    - end (int): End position of the sequence.

    Returns:
    - sequence (str): Fetched sequence in FASTA format.
    """
    # Provide your email address to NCBI
    #Entrez.email = email

    try:
        # Fetch the sequence from NCBI for the specified accession number and range
        print(f"fetching sequence for: {accession_number}")
        handle = Entrez.efetch(db="nuccore", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
        record = SeqIO.read(handle, "fasta")
        handle.close()

        sequence = str(record.seq)
        return sequence

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def fetch_query(query_fasta, start, end):
    """
    Fetch a specified range from a query sequence in a FASTA file.

    Parameters:
    - query_fasta (str): Path to the FASTA file containing the query sequence.
    - start (int): Start position of the desired range.
    - end (int): End position of the desired range.

    Returns:
    - sequence (str): Fetched sequence in the specified range.
    """
    try:
        # Read the query sequence from the FASTA file
        with open(query_fasta, "r") as file:
            records = list(SeqIO.parse(file, "fasta"))

        if not records:
            raise ValueError("No sequences found in the FASTA file.")

        # Assume there is only one sequence in the file
        query_sequence = records[0].seq

        # Fetch the specified range from the query sequence
        sequence_range = query_sequence[start - 1:end]

        return str(sequence_range)

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def merge_hit_seqs(hseqs_fasta, merge_type):
    fout = f"{hseqs_fasta.rsplit('.')[0]}_merged_by_{merge_type}.fasta"
    arguments = [hseqs_fasta, fout, merge_type]

    try:
        subprocess.run(['python', 'merge_webblast.py']+arguments, check=True)
        print(f"merged alignment created!: {fout}")
    except subprocess.CalledProcessError as e:
        print(f"Error merging: {e}")

def parse_blast(tree_root, query_fasta, merge_type, outfile):
    """
    Creates a nucleotide alignment of blast hits by positioning sequences where they hit the query in the blast search
    """
    blast_type = tree_root.find(".//BlastOutput_program").text

    if "TSA" in tree_root.find(".//Hit_def").text:
        db = "TSA"
    else:
        db = tree_root.find(".//BlastOutput_db").text.split(":")[0]


    blast_hits = []
    for hit_element in tree_root.findall(".//Hit"):
        blast_hit = BlastHit(hit_element, blast_type)
        blast_hits.append(blast_hit)

    # Write FASTA sequences to the output alignment file
    sequences = []
    query_length = get_sequence_length(query_fasta)
    for blast_hit in blast_hits:
        fasta_sequences = blast_hit.get_fasta_sequences(query_length)
        sequences+=fasta_sequences
            
    SeqIO.write(sequences, outfile,"fasta")

    # Merge hit sequences
    merge_hit_seqs(outfile, merge_type)







def main(xml_file, query_fasta, merge_type, outfile):
    # Parse XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()


    # Process blast hits
    parse_blast(root, query_fasta, merge_type, outfile)

    # Merge hit sequences




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BLAST XML results and create nucleotide alignments.")
    parser.add_argument("xml_file", help="Path to the BLAST XML file")
    parser.add_argument("query_fasta", help="Path to the query FASTA file")
    parser.add_argument("merge_type",choices=["id", "species"], help="Method for merging hit sequences.")
    parser.add_argument("outfile", help="Path to the output FASTA file")
    parser.add_argument("email", help="email address to use NCBI Entrez") 
    args = parser.parse_args()
    Entrez.email = args.email
    main(args.xml_file, args.query_fasta, args.merge_type, args.outfile)

