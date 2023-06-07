import openai
import os
import requests
import json
from Bio import Entrez
import requests
import xml.etree.ElementTree as ET
import html
import unicodedata
import re
import time

CHATGPT_MODEL = "gpt-3.5-turbo"
CHATGPT_URL =  "https://api.openai.com/v1/chat/completions"
# API endpoints
pubmed_api = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&id="
crossref_api = "https://api.crossref.org/works/"
opencitations_api = "https://opencitations.net/index/coci/api/v1/citations/"

# メールアドレスを設定してください（Entrezの利用条件）
Entrez.email = os.environ.get("EMAIL_ADDRESS")



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
    

    #PubmedID 1
#def get_pubmed_id(title):
#        # title が "Not title found"でない場合のみ、PubMed IDを取得する
#    if title == "No title found":
#        return "No title found"
#    
#    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
#    search_url = base_url + "esearch.fcgi"
#    params = {
#        'db': 'pubmed',
#        'term': f"""\"{title}\"[Title:~1]""",
#        'retmode': 'xml'
#    }
#    try:
#        response = requests.get(search_url, params=params)
#        root = ET.fromstring(response.content)
#        ids = [id_elem.text for id_elem in root.iter('Id')]
#        if ids:
#            return ids[0]
#        else:
#            return "Not found"
#    
#    except Exception as e:
#        print(f"Error retrieving PubMed ID for \"{title}\": {e}")
#        return "Error"

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
    
def normalize_title(title):
    # Unescape HTML markup
    title = html.unescape(title)
    title = re.sub(r'<[^>]*?>', '', title)
    title = title.replace('[','').replace(']','')
    # Normalize unicode characters
    title = unicodedata.normalize('NFKD', title)
    return title
    

def main():
# 同じディレクトリにあるpapers.txtファイルから論文を読み込む
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, 'papers.txt')
    pubmed_ids = []

# 各論文からタイトルを抽出し、そのタイトルを使ってPubMed IDを取得する

    with open(file_path, 'r', encoding = 'utf-8') as f:
        papers = f.read().split('\n')
        for paper in papers:
            #print(paper)
            if paper:  # 空行を無視する
                title = extract_paper_title_from_text(paper)
                title = normalize_title(title)
                pubmed_id = get_pubmed_id(title)
                pubmed_ids.append(pubmed_id)

# paper, pubmed_id, citation_countの結果をファイルに書き出す
    out_file_path = os.path.join(current_dir, 'pubmed_ids_gpt4.txt')
    idx = 1
    with open(out_file_path, 'w', encoding = 'utf-8') as f:
        for paper, pubmed_id in zip(papers, pubmed_ids):
            f.write(f"{idx}:: {pubmed_id}\n")
            idx = idx+1

if __name__ == "__main__":
    main()