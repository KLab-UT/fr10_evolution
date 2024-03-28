import argparse
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import re
import csv

class BlastHit:
    """
    Class to represent a sequence with hits from NCBI blast results
    A BlastHit can represent multiple alignments from a certain sequence to the query sequence, i.e., hsps or high scoring pairs
    """

    def __init__(self, hit_element):
        self.accession = hit_element.find("Hit_accession").text
        self.definition = hit_element.find("Hit_def").text
        self.length = int(hit_element.find("Hit_len").text)
        self.hsps = hit_element.findall("Hit_hsps/Hsp")
        self.avg_identity = self.get_hsps_avg_identity()


    def get_hsps_avg_identity(self):
        identites = [int(hsp.find("Hsp_identity").text) for hsp in self.hsps]
        avg = sum(identites) / len(identites)
        return avg

def parse_blast_xml_hits(blast_result_xml):
    """
    return a list of all blast hit elements in xml results file
    """

    # Parse XML file
    tree = ET.parse(blast_result_xml)
    root = tree.getroot()

    blast_type = root.find(".//BlastOutput_program").text

    # Parse out blast hits
    blast_hits = []
    for hit_element in root.findall(".//Hit"):
        blast_hit = BlastHit(hit_element)
        blast_hits.append(blast_hit)

    print(f"total sequences hit: {len(blast_hits)}")
    return blast_hits

def log_annotations(blast_hits, annotations_csv):
    """
    Write the accession ids and annotations for a list of BlastHits
    """

    # Store accession ids and annotations from BlastHits
    id_annotations = []
    for hit in blast_hits:
        id = hit.accession
        avg_identity = hit.avg_identity

        # Clean annotations 
        annotation = re.sub(r'\[.*?\]', '', hit.definition) # remove species
        annotation = re.sub(r'>.*', '', annotation)  # Remove anything following ">"
        annotation = re.sub(r'isoform.*', '', annotation)  # Remove anything following "isoform"
        annotation = annotation.replace(',','') # Remove commas
        annotation = annotation.replace("PREDICTED: ", '')  # Remove "PREDICTED: "
        annotation = annotation.replace("-like", "") # Remove '-like'
        annotation = annotation.replace("putative ", "") # remove 'putative'
        annotation = annotation.replace("probable ","") # remove 'probable


        # Specific changes, adjust as needed
        if "A-I" in annotation or "apolipoprotein IA" in annotation:
            annotation = "apolipoprotein A-I"

        if "apolipoprotein A2" in annotation or "Apolipoprotein A-II" in annotation or "apolipoprotein A-II" in annotation or "APOA2" in annotation or "apolipoprotein a-ii" in annotation or "apoD" in annotation:
            annotation = "apolipoprotein A-II"

        if "APOA4" in annotation or "apolipoprotein A-IV" in annotation or "apolipoprotein A4" in annotation:
            annotation = "apolipoprotein A-IV"

        if "apolipoprotein C-I" in annotation or "apolipoprotein C Ib" in annotation or "APOC1" in annotation:
            annotation = "apolipoprotein C-I"

        if "apolipoprotein E" in annotation:
            annotation = "apolipoprotein E"

        if "hypothetical protein" in annotation or "uncharacterized protein" in annotation or "unnamed protein" in annotation:
            annotation = "hypothetical protein/uncharacterized protein"
        
        if "uveal autoantigen" in annotation:
            annotation = "uveal autoantigen" 
        
        if "zinc finger protein 853" in annotation:
            annotation = "zinc finger protein 853"

        if "LS-12" in annotation or "type-IV" in annotation or "type-4" in annotation or "type IV" in annotation or "Type-4" in annotation or "AFP4" in annotation or "type 4" in annotation:
            annotation = "LS-12 type-IV ice structuring protein"

        if "epidermal growth factor receptor substrate 15" in annotation:
            annotation = "epidermal growth factor receptor substrate 15"

        if "FYVE" in annotation:
            annotation = "FYVE domain-containing protein"


        # GROUP ALL APOS TOGETHER
        # if "apo" in annotation:
        #     annotation = "Apo A-I/A-II/A-IV/C-I/E"

        

        #print(f"Annotation for {id}\n{annotation}")
        id_ann = (id, annotation, avg_identity)
        id_annotations.append(id_ann)

    # Write annotations
    with open(annotations_csv, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        header = ("accession", "annotation", "identity")
        csv_writer.writerow(header)
        for id_ann in id_annotations:
            csv_writer.writerow(id_ann)




def main():
    parser = argparse.ArgumentParser(description='Process XML blast results and log annotations.')
    parser.add_argument('blast_result_xml', type=str, help='Path to the XML file containing blast results')
    parser.add_argument('annotations_csv', type=str, help='Path to the CSV file to store annotations')

    args = parser.parse_args()

    blast_hits = parse_blast_xml_hits(args.blast_result_xml)
    log_annotations(blast_hits, args.annotations_csv)

if __name__ == "__main__":
    main()
