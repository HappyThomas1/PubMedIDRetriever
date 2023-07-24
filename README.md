# PubMedIDRetriever
This is a suite of scripts designed to retrieve PubMed IDs or citations.

## Retrieve_PMID_and_Citations.py
With this script, you can extract PubMed IDs and their corresponding citation counts from a reference list, regardless of its format. The script utilizes ChatGPT to identify the title of each reference, hence there's no need for a specific format, as long as the ChatGPT API can identify the titles. 

### Requirements

The script is written in Python 3 and requires the following Python packages:

- openai
- requests
- Bio

```
pip install requests biopython openai
```

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


## RIS_to_RIS.py

This Python program helps you to extract literature references from an input Research Information Systems (RIS) file, fetch the related data from PubMed using the Digital Object Identifier (DOI), and store the results in an output RIS file.

### Requirements

The script is written in Python 3 and requires the following Python packages:

- requests
- Bio
- rispy
- bs4
- pandas

Before running the script, please make sure these packages are installed. You can install them using pip:

```
pip install requests biopython rispy beautifulsoup4 pandas
```

You need to set your email address as an environment variable for the Entrez module. Use the variable name `EMAIL_ADDRESS` for this purpose.

### Usage

Just run the script in Python. The script will open a file dialog where you can choose the input RIS file. The script will then extract all DOIs from the RIS file, fetch the related data from PubMed, and store the results in an output RIS file in a new directory named "output_ris_[input file name]" in the same directory as the script.

```
python RIS_to_RIS.py
```

### Workflow

The script follows these steps:

1. Opens a dialog to select the input RIS file.
2. Creates an output directory in the same location as the script.
3. Checks if there is a previous run that already fetched some DOIs and their data. If so, it loads these DOIs to skip fetching their data again.
4. Parses the input RIS file and extracts all DOIs.
5. For each DOI in the input file, the script fetches the citing articles and adds their DOIs to the list.
6. The script then fetches the PubMed data for each DOI in the list (excluding the ones already fetched in a previous run) and writes the data to the output RIS file.


The output is a RIS file named "combined.ris" stored in a directory named "output_ris_[input file name]" in the same location as the script. The RIS file contains the PubMed data for all DOIs extracted from the input file and their citing articles.

The script also stores all DOIs it has processed in a CSV file named "all_dois.csv" in the same directory as the output file. This is to avoid fetching the data for the same DOI again if the script is run multiple times with the same input file.



## BibTex_reformat.py

This Python script is designed to convert BibTex files into a specified format. It takes a BibTex file as an input, parses it, and produces a formatted text file. It has a built-in GUI that you can use to select the file and an optional reverse order function.

### Requirements

This script is written in Python 3 and requires the following Python packages:

- bibtexparser
- tkinter
- pylatexenc

```
pip install bibtexparser tkinter pylatexenc
```

### Usage

To use this script, simply run it using Python. It will open a window with a "Open BibTex File" button that you can use to select your BibTex file. If you want the resulting formatted bibliography to be in reverse order, you can check the 'Reverse Order' checkbox.

```
python BibTex_to_Format.py
```

### Workflow

1. Run the script to open the GUI.
2. Click the "Open BibTex File" button and select your BibTex file from the dialog box.
3. (Optional) If you want the resulting bibliography to be in reverse order, check the 'Reverse Order' box.
4. The script will parse the BibTex file, format the entries according to the specified style, and write them to a file named 'formatted_bibliography.txt' in the same directory as the input file.

This script is also capable of handling certain special LaTeX markups in the BibTex file, such as Greek alphabet symbols and superscripts.



