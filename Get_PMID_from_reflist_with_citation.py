import openai
import os
import requests
import json
from Bio import Entrez
import time

CHATGPT_MODEL = "gpt-3.5-turbo"
CHATGPT_URL =  "https://api.openai.com/v1/chat/completions"
# API endpoints
pubmed_api = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
crossref_api = "https://api.crossref.org/works/"
opencitations_api = "https://opencitations.net/index/coci/api/v1/citations/"

# メールアドレスを設定してください（Entrezの利用条件）
Entrez.email = os.environ.get("EMAIL_ADDRESS")


def get_paper_title_and_doi(pmid):
    # pmidが数字出ない場合は、"Error"を返す
    if not pmid.isdigit():
        return "Error"
    
    # Send a request to the PubMed API
    response = requests.get(pubmed_api + str(pmid))
    response.raise_for_status()  # Ensure we got a valid response
    paper_data = response.json()

    # Extract the paper title and DOI
    #paper_title = paper_data['result'][str(pmid)]['title']
    paper_doi = paper_data['result'][str(pmid)]['elocationid'].replace('doi: ', '')

    return  paper_doi

def get_citation_count(doi):
    # doiが"Error"の場合は、"Error"を返す
    if doi == "Error":
        return "Error"
    # Send a request to the OpenCitations API
    response = requests.get(opencitations_api + doi)
    response.raise_for_status()  # Ensure we got a valid response
    citation_data = response.json()

    # Extract the citation count
    citation_count = len(citation_data)  # Each item in the list is a citation

    return citation_count

def extract_paper_title_from_text(text):
    api_url = CHATGPT_URL
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {openai.api_key}"
    }
    data = {
        "model": CHATGPT_MODEL,
        "messages":  [{"role": "user", "content": "I have the following reference:\n\n" + text + "\n\nThe title of the paper is:"}],
        "max_tokens": 100,
        "n": 1,
        "stop": None,
        "temperature": 0,
    }
    print(text)

    response = requests.post(api_url, headers=headers, data=json.dumps(data))
    print(response)
    result = response.json()
    
    if 'choices' in result and len(result['choices']) > 0:
        title = result["choices"][0]["message"]["content"].strip()
        print(title)
        return title
    else:
        return "No title found"

def get_pubmed_id(title):
    # title が "Not title found"でない場合のみ、PubMed IDを取得する
    if title == "No title found":
        return "No title found"
    # title が "Error"でない場合のみ、PubMed IDを取得する
    max_retries = 3
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="pubmed", term=f"""\"{title}\"[Title:~1]""", retmax=1)
            record = Entrez.read(handle)
            handle.close()
            #time.spleep(1)

            if record["IdList"]:
                return record["IdList"][0]
            else:
                return "Not found"
       
        except Exception as e:
            print(f"Error retrieving Title for \"{title}\" on attempt {attempt+1}: {e}")
            time.sleep(1)
    
    print(f"Failed to retrieve Title for \"{title}\" after {max_retries} attempts.")
    return "Error"
    

def main():
# 同じディレクトリにあるpapers.txtファイルから論文を読み込む
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, 'papers.txt')
    pubmed_ids = []
    doi_list = []
    citation_count_list = []

# 各論文からタイトルを抽出し、そのタイトルを使ってPubMed IDを取得する

    with open(file_path, 'r',) as f:
        papers = f.read().split('\n')
        for paper in papers:
            #print(paper)
            if paper:  # 空行を無視する
                title = extract_paper_title_from_text(paper)
                pubmed_id = get_pubmed_id(title)
                pubmed_ids.append(pubmed_id)
                try:
                    doi = get_paper_title_and_doi(pubmed_id)
                    citation_count = get_citation_count(doi)
                
                    doi_list.append(doi)
                    citation_count_list.append(citation_count)
                
                except Exception as e:
                    print(f"An error occurred: {e}")

# paper, pubmed_id, citation_countの結果をファイルに書き出す
    out_file_path = os.path.join(current_dir, 'pubmed_ids.txt')
    with open(out_file_path, 'w') as f:
        for paper, pubmed_id, citation_count in zip(papers, pubmed_ids, citation_count_list):
            f.write(f"{paper}: {pubmed_id}, {citation_count}\n")

if __name__ == "__main__":
    main()