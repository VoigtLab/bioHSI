# Biomolecular network
This repository contains the code needed to extract metabolite and metabolic 
reaction information from the [BKMS](https://bkms.brenda-enzymes.org/), [Rhea](https://www.rhea-db.org/), and [BioCyc](https://biocyc.org/) databases

- `00_parse_reaction_dbs`: code to extract metabolites, reactions, and reaction information from the reaction databases.
  - Data: database flatfiles are provided at `bkms/`,`metacyc/`, and `rhea/` 

- `01_construct_network_of_biosynthesis`: code to combine information from the reaction databases and generate a network of biochemistry

- `02_parse_organism_dbs`: define host-specific metabolism lists
