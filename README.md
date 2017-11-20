# scvep
Super Cereal VEP

Annotation-free variant effect prediction for, well, unannotated reference genomes/transcriptomes. 

Requires ktio library: pip install ktio (--target /path/to/scvep *IFF* no root/superuser access/own python environment on system) and a local blast+ installation, i.e. blastx in path

===Installation===
git clone https://github.com/krasileva-group/scvep.git
cd scvep
pip install -e .


1. Build reference database for target organism/family etc., according to the example in prepare_db.sh
2. Run with python scvep.py <snps.vcf> <scaffolds.fasta> </path/to/reference_db>. scvep will generate two tables, a) unfiltered blast output and b) predictions on filtered blast hits.
