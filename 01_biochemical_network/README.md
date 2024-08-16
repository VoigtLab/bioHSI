# Biomolecular network
This repository contains the code needed to extract metabolite and metabolic 
reaction information from the [BKMS](https://bkms.brenda-enzymes.org/), [Rhea](https://www.rhea-db.org/), and [BioCyc](https://biocyc.org/) databases

- `01_parse_reaction_dbs`: code to extract metabolites, reactions, and reaction information from the reaction databases.
  - `00_data`: Raw input files and processed output files 
      - Enzymatic reaction flatfile directories: `00_data/raw/bkms/`, `00_data/raw/metacyc/`, `00_data/raw/rhea/`
  - `01_parse_reaction_dbs`: Code to clean, standardize, and retrieve information for the reaction databases
      - `query_uniprot_for_metacyc.py`: Retrieves enzyme sequences for MetaCyc reactions from uniprot 
      - `scrape_brenda.py`: 


- `02_construct_network_of_biosynthesis`: code to combine information from the reaction databases and generate a network of biochemistry

- `03_parse_organism_dbs`: define host-specific metabolism lists
