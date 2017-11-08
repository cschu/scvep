# scvep
Super Cereal VEP

Annotation-free variant effect prediction for, well, unannotated reference genomes/transcriptomes. 

Requires ktio library: pip install ktio and a local blast+ installation, i.e. blastx in path

1. Build reference database for target organism/family etc., according to the example in prepare_db.sh
2. Run with python scvep.py <snps.vcf> <scaffolds.fasta> </path/to/reference_db>. scvep will generate two tables, a) unfiltered blast output and b) predictions on filtered blast hits.
