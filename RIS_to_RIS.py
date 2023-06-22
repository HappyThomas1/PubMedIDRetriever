import os
import requests
from Bio import Entrez
import rispy
from bs4 import BeautifulSoup
import tkinter as tk
from tkinter import filedialog
import time
import pandas as pd


Entrez.email = os.getenv("EMAIL_ADDRESS")  # Please set your email address here.
RETRY_LIMIT = 3
RETRY_DELAY = 1  # エラー発生後に次のリクエストを送信するまでの待ち時間(秒)

def parse_ris(file_path):
    with open(file_path, 'r', encoding = 'utf-8') as bibliography_file:
        entries = list(rispy.load(bibliography_file))
    return entries



def get_pubmed_id(doi):
    for _ in range(RETRY_LIMIT):
        try:
            handle = Entrez.esearch(db='pubmed', term=doi)
            record = Entrez.read(handle)
            handle.close()
            if record["Count"] == "0":
                print("No results found for DOI.")
                return None
            pubmed_id = record["IdList"][0]
            print(pubmed_id)
            return pubmed_id
        except Exception as e:
            print(f"Error occurred on {doi}: {str(e)}, retrying...")
            time.sleep(RETRY_DELAY)
    return None

def fetch_from_pubmed(pubmed_id):
    for _ in range(RETRY_LIMIT):
        try:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pubmed_id}&retmode=xml"
            response = requests.get(url, verify=True)
            soup = BeautifulSoup(response.content, "xml")
            return soup
        except Exception as e:
            print(f"Error occurred on {pubmed_id}: {str(e)}, retrying...")
            time.sleep(RETRY_DELAY)
    return None

def fetch_references(doi):
    for _ in range(RETRY_LIMIT):
        try:
            url = f"https://api.crossref.org/works/{doi}"
            response = requests.get(url)
            work = response.json()['message']
            references = work.get('reference', [])            
            return references or []
        except Exception as e:
            print(f"Error occurred on {doi}: {str(e)}, retrying...")
            time.sleep(RETRY_DELAY)
    return []

def fetch_pubmed_data(doi):
    ris_data = None
    for _ in range(RETRY_LIMIT):
        try:
            pubmed_id = get_pubmed_id(doi)
            soup = fetch_from_pubmed(pubmed_id)
            if soup is not None:
                ris_data = convert_to_ris(soup, pubmed_id)
            
            return ris_data
        except Exception as e:
            print(f"Error occurred: {str(e)}, retrying...")
            time.sleep(RETRY_DELAY)
    return None

def convert_to_ris(soup, pubmed_id):
    ris = "TY  - JOUR\n"
    try:
        title = soup.find("ArticleTitle").text
        ris += f"T1  - {title}\n"
        ris += f"TI  - {title}\n"
    except AttributeError:
        print(f"No ArticleTitle for PubMed ID: {pubmed_id}")
        ris += f"T1  - ""\n"
        ris += f"TI  - ""\n"

    ris += f"AN  - {pubmed_id}\n"

    
    try:
        authors = soup.find_all("Author")
        for author in authors:
            firstname = author.ForeName.text
            lastname = author.LastName.text
            ris += f"AU  - {lastname}, {firstname}\n"
    except:
        ris += f"AU  - ""\n"
        
    try:
        year = soup.find("PubDate").Year.text
        ris += f"PY  - {year}\n"
    except:
        ris += f"PY  - ""\n"

    try:
        journal = soup.find("Title").text
        ris += f"JO  - {journal}\n"
    except:
        ris += f"JO  - ""\n"

    try:
        abstract = soup.find("Abstract").text
        ris += f"AB  - {abstract}\n"
    except:
        ris += f"AB  - ""\n"
    
    try:
        doi = soup.find("ELocationID").text
        ris += f"DO  - {doi}\n"
    except:
        ris += f"DO  - ""\n"

    try:
        vol = soup.find("Volume").text
        ris += f"VL  - {vol}\n"
    except:
        ris += f"VL  - ""\n"


    ris += f"UR  - https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/\n"
    ris += f"M1  - {pubmed_id}\n"
    ris += "ER  - \n"
    return ris

def create_output_directory(input_file_path):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_name = os.path.basename(input_file_path)
    file_name = os.path.splitext(file_name)[0]
    output_dir = os.path.join(current_dir, "output_ris_"+ file_name )
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    return output_dir

def select_ris_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(filetypes=[('RIS files', '*.ris')])
    root.destroy()
    return file_path

def get_existing_dois(file_path):
    """Read the existing DOIs from the given RIS file."""
    if not os.path.exists(file_path):
        return set()

    ris_entries = parse_ris(file_path)
    existing_dois = {entry.get('doi') for entry in ris_entries if entry.get('doi')}
    return existing_dois

# In your main function
def main():
    file_path = select_ris_file()
    max_depth = 1
    idx = 0
    if file_path:
        # Create output directory and get output file path
        output_dir = create_output_directory(file_path)
        output_file_path = os.path.join(output_dir, "combined.ris")

        # Load DOIs if they exist
        all_dois_file_path = os.path.join(output_dir, "all_dois.csv")
                    # Get existing DOIs in the output file
        existing_dois = get_existing_dois(output_file_path)

        if os.path.exists(all_dois_file_path):
            all_dois = set(pd.read_csv(all_dois_file_path, header=None)[0])
        else:
            # Parse the original RIS file
            ris_entries = parse_ris(file_path)
            all_dois = {entry.get('doi') for entry in ris_entries if entry.get('doi')}

            # Fetch DOIs and PubMed data
            for original_doi in all_dois.copy():
                depth = 0
                doi_list = [original_doi]
                while doi_list and depth < max_depth:
                    next_doi_list = []
                    for doi in doi_list:
                        references = fetch_references(doi)
                        for reference in references:
                            ref_doi = reference.get('DOI')
                            # Skip if the DOI is already in the output file
                            if ref_doi in existing_dois:
                                continue
                            if ref_doi and ref_doi not in all_dois:
                                print(f"""{idx}. DOI: {ref_doi}""")
                                idx += 1
                                next_doi_list.append(ref_doi)
                                all_dois.add(ref_doi)
                    doi_list = next_doi_list
                    depth += 1

            # Save all_dois to a file
            pd.DataFrame(list(all_dois)).to_csv(all_dois_file_path, index=False, header=False)

        print(all_dois)

        # Fetch and write PubMed data for new DOIs
        for doi in all_dois - existing_dois:
            ris_data = fetch_pubmed_data(doi)
            if ris_data:
                with open(output_file_path, 'a', encoding='utf-8') as f:
                    f.write(ris_data)

if __name__ == "__main__":
    main()
