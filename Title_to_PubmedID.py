import os
from Bio import Entrez
import csv
import html
import unicodedata
import re
import time


# Your email here.
# Entrez requires you to specify your email address with each request.
Entrez.email = os.environ.get("EMAIL_ADDRESS")

#def get_pubmed_id_from_title(title):
#    handle = Entrez.esearch(db="pubmed", term=title)
#    record = Entrez.read(handle)
#    return record["IdList"][0]  # Returns the first matching PubMed ID

def normalize_title(title):
    # Unescape HTML markup
    title = html.unescape(title)
    title = re.sub(r'<[^>]*?>', '', title)
    title = title.replace('[','').replace(']','')
    # Normalize unicode characters
    title = unicodedata.normalize('NFKD', title)
    return title

def get_pubmed_id_from_title(title):
    print(title)
    max_retries = 3
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="pubmed", term=f"""\"{title}\"[Title:~1]""", retmax=1)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(1)
            
            if record["IdList"]:
                return record["IdList"][0]
            else:
                return "Not found"
            
        except Exception as e:
            print(f"Error retrieving Title for \"{title}\" on attempt {attempt+1}: {e}")
            time.sleep(1)
    
    print(f"Failed to retrieve Title for \"{title}\" after {max_retries} attempts.")
    return "Error"

def read_and_write(input_file, output_file):
    with open(input_file, 'r', encoding = 'utf-8') as in_file, open(output_file, 'w',encoding = 'utf-8') as out_file:
        for line in in_file:
            row = line.strip().split('::')
            title = row[2].strip()
            #title = row[1].split(',')[1:]
            #title = ''.join(title)
            #title = title.strip()
            title = normalize_title(title)
            pubmed_id = get_pubmed_id_from_title(title)
            print(pubmed_id)
            out_file.write(f"{row[0]}:: {pubmed_id}:: {title}\n")


def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir,'pubmed_ids_title.txt') # Your input file here.
    output_file =os.path.join(current_dir,'pubmed_ids_newsearch_condition.txt') # Your output file here.

    read_and_write(input_file, output_file)


if __name__ == '__main__':
    main()