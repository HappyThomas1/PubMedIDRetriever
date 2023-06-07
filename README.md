# GetPubMedID
Collection of scripts that can obtain pubmed ID or citation

## Get_PMID_from_reflist_with_citation.py
If you have a reference list with some format, the program obtain pubmed ID along with the number of citation index.
As the program use ChatGPT to extract the title of the paper to request, the reference list does not need formatted as long as 
ChatGPT API recognize which is the title of the reference. 

Requirement:
OpenAI API key should be set as an environmental variable OPENAI_API_KEY
Your email address set as an environmental variable EMAIL_ADDRESS. Alternatively, you may directly edit the code as follows
```
Entrez.email = 'Your email address'
```
## Get_PMID_from_reflist.py
Almost same as Get_pubmedID_from_ref_list_with_citation.py, but only return the PubMedID list instead of citation index.

## PMID_to_Title.py
If you give a list of pubmedID, the script returns the list of the titles of the reference.

input file should be named as "pubmed_ids.txt"

Format of pubmed_ids.txt
```
1:: PMID_1
2:: PMID_2
```

## Title_to_PMID.py
If you give a list of titles of papers, the script returns the list of pubmedIDs of the reference.

input file should be named as "pubmed_title.txt"

Format of pubmed_title.txt
```
1:: title_1
2:: title_2
```


## PMID_to_RIS.py
The script assumes word file using EndNote plugin and you have references brancketed like '{PMID}'. Drag and Drop the word file, then
the script extract all the possible PMIDs and return the RIS format for reference.

