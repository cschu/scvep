import sys
import os
import argparse

from scvep import *

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('vcf', type=str)
    ap.add_argument('reference', type=str)
    ap.add_argument('blastdb', type=str)
    ap.add_argument('--flanksize', '-f', type=int, default=100)
    ap.add_argument('--prefix', '-p', type=str, default='scvep')
    args = ap.parse_args()
    assert os.path.exists(args.vcf)
    assert os.path.exists(args.reference)
    assert os.path.exists(args.blastdb)
    assert args.flanksize > 50

    with open(args.vcf) as vcf_in:
        print("Reading SNPs...")
        snps = readSNPs(vcf_in)
    with open(args.reference) as seq_in:
        print("Processing SNPs...")
        processSNPs(seq_in, args)    
