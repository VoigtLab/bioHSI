# Parse Reaction Databases
This directory contains the code needed to parse and clean enzymatic reaction data (chemical transformation and enzyme sequences) from reaction databases.

Data is retrieved and saved to `../00_data` by default.

### Description of code in this directory
The processing scripts are largely organized by the database they are used for.
For KEGG, Brenda, and MetaCyc, the BKMS dataset 

##### MetaCyc and Rhea
The notebook `parse_rhea_and_metacyc.ipynb` contains code to process the data retrieved from the MetaCyc and Rhea databases, convert the reactions into SMILES, clean the reaction SMILES, and remove co-factors as needed. 

`query_uniprot_for_metacyc.py` and `query_uniprot_for_rhea.py` contain code to submit API requests to the Uniprot database to retrieve enzyme sequences for provided protein IDs. 

The notebook `retrieve_sequences_for_metacyc_and_rhea.ipynb` contains code to collate reaction and protein sequence information

##### Brenda

`scrape_brenda.py` contains code to extract Reaction and protein accession number information from the Brenda database. 

The outputs are used in `query_uniprot_for_brenda.py` to retrieve enzyme sequences for Brenda reactions.

##### KEGG

`query_kegg_for_bkms.py` contains code to scrape KEGG based on the reaction data in the BKMS database. Retrieves mapping from KEGG reactions to sequence accession enzyme sequences

##### Compilation in BKMS
The data retrieved for MetaCyc, Brenda, and KEGG are compiled with teh BKMS database in `retrieve_sequences_for_bkms.ipynb`