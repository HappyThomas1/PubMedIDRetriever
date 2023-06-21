import os
import sys
import requests
from bs4 import BeautifulSoup
import re
from docx import Document
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinterdnd2 import TkinterDnD
from requests.exceptions import RequestException


def extract_arabic_numerals(docx_path):
    doc = Document(docx_path)
    pattern = r'\{(\d+)\}'
    arabic_numerals = []

    for paragraph in doc.paragraphs:
        matches = re.findall(pattern, paragraph.text)
        arabic_numerals.extend(matches)

    return arabic_numerals

def fetch_pubmed_article(pubmed_id):
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pubmed_id}&retmode=xml"
        response = requests.get(url)
        
        # Raise an error except when requests.get was successful
        if response.status_code != 200:
            raise ValueError(f"PubMed ID: {pubmed_id} not found.")
        
        soup = BeautifulSoup(response.content, "xml")
        
        # Raise an error when the value was empty or not found
        if not soup or not soup.find('PubmedArticle'):
            raise ValueError(f"PubMed ID: {pubmed_id} not found or response not as expected.")
        
        return soup

    except RequestException as e:
        print(f"Error occurred while fetching the data: {e}")
        return None

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None



# 
def convert_to_ris(soup, pubmed_id):
    
    ris = "TY  - JOUR\n"
    title = soup.find("ArticleTitle").text
    ris += f"T1  - {title}\n"
    ris += f"AN  - {pubmed_id}\n"

    authors = soup.find_all("Author")
    for author in authors:
        lastname = author.LastName.text
        try:
            firstname = author.ForeName.text
        except:
            firstname = ""
        ris += f"AU  - {lastname}, {firstname}\n"

    year = soup.find("PubDate").Year.text
    ris += f"PY  - {year}\n"

    try:
        journal = soup.find("Title").text
        ris += f"JO  - {journal}\n"
    except:
        pass

    try:
        abstract = soup.find("Abstract").text
        ris += f"AB  - {abstract}\n"
    except:
        pass
    
    try:
        doi = soup.find("ELocationID").text
        ris += f"DO  - {doi}\n"
    except:
        pass

    try:
        vol = soup.find("Volume").text
        ris += f"VL  - {vol}\n"
    except:
        pass


    ris += f"UR  - https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/\n"
    ris += f"M1  - {pubmed_id}\n"
    ris += "ER  - \n"
    return ris


def on_drop(event):
    filename = event.data.strip("{}")
    if filename.lower().endswith(".docx"):
        #print(filename)
        pubmed_ids = extract_arabic_numerals(filename)
        # Get the directory of the script
        script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        if not os.path.exists(os.path.join(script_dir, "Outputs")):
            os.mkdir(os.path.join(script_dir, "Outputs"))
        
        script_dir = os.path.join(script_dir, "Outputs")
        print(script_dir)
        # Create the output file path in the same directory as the script
        output_file = os.path.split(filename)[1]
        output_file = os.path.splitext (output_file) [0]
        output_path = os.path.join(script_dir, output_file+ "_output.ris")
        print(output_path)
        with open(output_path, "w",encoding= 'utf-8') as f:
            for pubmed_id in pubmed_ids:
                article_soup = fetch_pubmed_article(pubmed_id)
                
                # If article_soup is None, skip this pubmed_id
                if article_soup is None:
                    print(f"Skipping PubMed ID: {pubmed_id} due to an error while fetching the data.")
                    continue

                ris = convert_to_ris(article_soup, pubmed_id)
                f.write(ris)
        print(f"PubMed IDs extracted and output.ris file generated.")
    else:
        print("Invalid file format. Please drop a .docx file.")


def main():
    global label
    # Create the main window using TkinterDnD
    root = TkinterDnD.Tk()
    root.geometry('400x200')
    root.title("PubMed ID Extractor")
    
    label = tk.Label(root, text="Drag and Drop Word Document", font=("Helvetica", 18))
    label.pack(padx=20, pady=20)

    tkdnd_path = os.path.join("C:", os.sep, "Program Files", "Python311", "tcl", "tkdnd2.8")
    root.tk.call("lappend", "auto_path", tkdnd_path)

    root.drop_target_register('*')
    root.dnd_bind('<<Drop>>', on_drop)

    root.mainloop()

if __name__ == "__main__":
    main()