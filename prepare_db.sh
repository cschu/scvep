#!/bin/bash

# NCBI Genbank query: "Poaceae"[Organism] NOT chloroplast[filter] NOT mitochondrion[filter] NOT plastid[filter] NOT plasmid[filter] 
# TODO: Add Entrez query here

source cdhit-4.6.1
sbatch -p tgac-medium -c 8 --mem 32GB <(echo '#!/bin/bash -e'$'\n'"cd-hit -i grass_proteins.fa -o grass_proteins.cd95.fa -c 0.95 -T 8 -M 32000")

source blast-2.6.0
makeblastdb -in grass_proteins.cd95.fa -dbtype prot -parse_seqids
