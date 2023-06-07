import os
from Bio import Entrez
import csv
import time

# Your email here.
# Entrez requires you to specify your email address with each request.
Entrez.email = os.environ.get("EMAIL_ADDRESS")

def get_title_from_pubmed_id(pubmed_id):
    max_retries = 3
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml",retmax=1)
            records = Entrez.read(handle)
            handle.close()
            # Add a delay of 1 second (or more if needed)
            time.sleep(1)
            return records['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
        
        except Exception as e:
            print(f"Error retrieving Title for \"{pubmed_id}\" on attempt {attempt+1}: {e}")
            time.sleep(1)
    
    print(f"Failed to retrieve Title for \"{pubmed_id}\" after {max_retries} attempts.")
    return "Error"

    

def read_and_write(input_file, output_file):
    with open(input_file, 'r', encoding = 'utf-8') as in_file, open(output_file, 'w', encoding = 'utf-8') as out_file:
        for line in in_file:
            row = line.strip().split('::')
            pubmed_id = row[1].strip()
            title = get_title_from_pubmed_id(pubmed_id)
            print(title)
            out_file.write(f"{row[0]}:: {pubmed_id}:: {title}\n")


def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir,'pubmed_ids.txt') # Your input file here.
    output_file =os.path.join(current_dir,'pubmed_ids_title.txt') # Your output file here.

    read_and_write(input_file, output_file)


if __name__ == '__main__':
    main()