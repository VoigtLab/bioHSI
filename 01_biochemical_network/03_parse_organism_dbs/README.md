Download appropriate data from biocyc website
1. Find BioCyc webiste or database specific to organism of interest (e.g. https://yeast.biocyc.org/?sid=biocyc16-3884907109 for yeast)
2. Navigate to Tools > SmartTables > Special SmartTables
3. Open the "All Reactions of X" table and make a personal, editable copy.
4. Choose Add Property > In-pathway and Substrates
5. Download Export > to Spreadsheet file... (with common names)
6. Go back to Special SmartTables
7. Get the SmartTable "All compounds of X"
8. Add the SMILES property Column
9. Export 
10. Use `extract_organism_metabolites.ipynb` to extract the metabolites in a format that can be used for the biosynthesis network search. 

A similar process can be followed from MetanetX for which `extract_organism_metabolites_metanetx.ipynb` can be used.