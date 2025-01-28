# Parse Reaction Databases
This directory contains the code needed to parse and clean enzymatic reaction data (chemical transformation and enzyme sequences) from reaction databases.

Data is retrieved and saved to `../00_data` by default.

### Description of code in this directory
The processing scripts are largely organized by the database they are used for.

##### MetaCyc and Rhea
The notebook `parse_rhea_and_metacyc.ipynb` contains code to process the data retrieved from the MetaCyc and Rhea databases, convert the reactions into SMILES, clean the reaction SMILES, and remove co-factors as needed. 

The `query_uniprot<...>.py` scripts contain code to submit API requests to the Uniprot database to retrieve enzyme sequences for provided protein IDs. 

The notebook `retrieve_sequences_for_metacyc_and_rhea.ipynb` contains code to collate reaction and protein sequence information
