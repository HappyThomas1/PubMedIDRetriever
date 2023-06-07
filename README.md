# PubMedIDRetriever
This is a suite of scripts designed to retrieve PubMed IDs or citations.

## Retrieve_PMID_and_Citations.py
With this script, you can extract PubMed IDs and their corresponding citation counts from a reference list, regardless of its format. The script utilizes ChatGPT to identify the title of each reference, hence there's no need for a specific format, as long as the ChatGPT API can identify the titles. 

### Requirements:
Your input file should be named 'papers.txt', with each line containing a single reference.

Set your OpenAI API key as an environmental variable named OPENAI_API_KEY.

Set your email address as an environmental variable named EMAIL_ADDRESS. Alternatively, you can directly modify the code as follows:
```
Entrez.email = 'Your email address'
```

## Retrieve_PMID.py
This script is similar to Retrieve_PMID_and_Citations.py, but it only returns a list of PubMed IDs, excluding the citation counts.

## PMID_to_Title.py
Provide a list of PubMed IDs and this script will return the corresponding titles for each reference.

Your input file should be named "pubmed_ids.txt".

The format for pubmed_ids.txt should be:
```
1:: PMID_1
2:: PMID_2
```

## Title_to_PMID.py
Supply a list of paper titles and this script will return the associated PubMed IDs.

Your input file should be named "pubmed_title.txt".

The format for pubmed_title.txt should be:
```
1:: title_1
2:: title_2
```

## PMID_to_RIS.py
This script is designed for Word files using the EndNote plugin, where references are bracketed like '{PMID}'. Simply drag and drop the Word file, and the script will extract all possible PubMed IDs and convert them into the RIS reference format.
