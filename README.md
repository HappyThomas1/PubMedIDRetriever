# GetPubMedID
Collection of scripts that can obtain pubmed ID or citation

## Get_pubmedID_from_ref_list_with_citation
If you have a reference list with some format, the program obtain pubmed ID along with the number of citation index.
As the program use ChatGPT to extract the title of the paper to request, the reference list does not need formatted as long as 
ChatGPT API recognize which is the title of the reference. 

Requirement:
OpenAI API key should be set as an environmental variable OPENAI_API_KEY
Your email address set as an environmental variable EMAIL_ADDRESS. Alternatively, you may directly edit the code as follows
```
Entrez.email = 'Your email address'
```
