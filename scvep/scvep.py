#!/usr/bin/env python
# coding=utf-8
import sys
import csv
import subprocess as sub
import os
import argparse
from collections import namedtuple

from ktio.ktio import readFasta

__author__ = 'Christian Schudoma (cschu)'
__copyright__ = 'Copyright 2017, Christian Schudoma, Earlham Institute'
__license__ = 'MIT'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma (cschu)'
__email__ = 'cschu1981@gmail.com'

# BLAST_DB = os.path.join(os.path.dirname(sys.argv[0]), 'grass_proteins.cd95.fa')
# BLAST_CMD = 'blastx -max_target_seqs 1 -db {} -outfmt "6 std sseqid qstart qend sstart send evalue bitscore qlen slen positive gaps ppos qframe staxids salltitles qseq sseq"'
BLAST_CMD = 'blastx -max_target_seqs 1 -db {} -outfmt "6 std sseqid qlen slen positive gaps ppos qframe salltitles qseq sseq"'
MARKER_LENGTH = 201

VariantPosition = namedtuple('VariantPosition', 'seqid pos ref alt'.split(' '))
BlastHSP = namedtuple('BlastHSP', 'query subject pident length mismatch gapopen qstart qend sstart send evalue bitscore sseqid qlen slen positive gaps ppos qframe salltitles qseq sseq'.split(' '))

codonWheel = { 'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
               'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
               'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
               'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
               'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
               'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
               'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
               'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
               'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
               'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
               'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
               'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
               'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
               'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
               'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
               'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F', }

def translate(codon):
    return codonWheel.get(codon.upper().replace('U', 'T'), 'X')

def reverseComplement(seq):
    """
    Returns the reverse complement of a nucleic acid sequence input.
    """
    compl= dict(zip('ACGTNRYWSMKBHDV', 'TGCANYRWSKMVDHB'))
    return ''.join([compl[base]
                    for base in seq.upper().replace('U', 'T')])[::-1]

def filterBlastHit(blast_hit, pos_threshold=50, qcov_threshold=0.6, region_length=MARKER_LENGTH):
    return int(blast_hit.gaps) == 0 and float(blast_hit.ppos) >= pos_threshold and (3 * float(blast_hit.length)) / region_length >= qcov_threshold

def predictEffect(seq, blast_hit, ref_allele, alt_allele, _out):
    def extractCodon(seq, strand='+'):
        """
        due to filtering for qcov > 50% (higher, better), we can
        use the length of the blast hit on the left flank to determine
        the codon position of the SNP - assuming that the blast hit
        represents the actual identity of the region
        """
        cpos = len(seq[0]) % 3
        if cpos == 0:
            """ |left flank - qstart| is divisible by 3, i.e. SNP starts a new codon """
            # print(seq[0], seq[1] + seq[2][:2], seq[2][2:], sep='_')
            return seq[1] + seq[2][:2], cpos
        elif cpos == 1:
            """ new codon starts at last position of left flank, SNP center """
            # print(seq[0][:-1], seq[0][-1] +  seq[1] + seq[2][0], seq[2][1:], sep='_')
            return seq[0][-1] + seq[1] + seq[2][0], cpos
        elif cpos == 2:
            """ SNP is at position 3 of codon """
            # print(seq[0][:-2], seq[0][-2:] + seq[1], seq[2], sep='_')
            return seq[0][-2:] + seq[1], cpos

    frame = int(blast_hit.qframe)
    qstart, qend = sorted(map(int, (blast_hit.qstart, blast_hit.qend)))
    # varseq = seq[0] + seq[1] + seq[2]

    if qend < qstart:
        """
         deal with negative strand hits, i.e. swap qstart <-> qend
        """
        qstart, qend = qend, qstart
        # varseq = reverseComplement(seq[2]) + reverseComplement(seq[1]) + reverseComplement(seq[0])

    """
     filtering assures that we only see regions here that are covered by blast hit
     hence, qstart is located in left flank and qend in right flank, i.e. SNP is covered as well
    """
    seq, strand = (seq[0][qstart-1:], seq[1], seq[2][:qend-len(seq[0])-1]), '+'
    if frame < 0:
        """
         now take reverse strand if negative frame was hit
        """
        seq, strand = tuple(reversed(list(map(reverseComplement, seq)))), '-'        
        ref_allele, alt_allele = map(reverseComplement, [ref_allele, alt_allele])

    # print(*seq, sep='_')

    assert ref_allele == seq[1]
    """ codon extraction s. above """
    codon, cpos = extractCodon(seq, strand)
    alt_codon = list(codon)
    alt_codon[cpos] = alt_allele
    alt_codon = ''.join(alt_codon)
    # blast_hit[20]
    out = [blast_hit.query, blast_hit.sseqid, cpos, frame, 'bh[20]', blast_hit.pident, abs(float(blast_hit.qend)-float(blast_hit.qstart))/float(blast_hit.qlen), blast_hit.evalue]
    mclass = 'nonsense' if translate(alt_codon) == '*' else ('missense' if translate(codon) != translate(alt_codon) else 'silent')
    out.extend([codon, translate(codon), alt_codon, translate(alt_codon), mclass])
    print(*out, sep='\t', file=_out, flush=True)
    # print(*out, sep='\t')


def readVariantPositions(_in):
    for row in csv.reader(_in, delimiter='\t'):
        if not row[0].startswith('#'):
            yield row

def readSNPs(_in):
    snps = dict()
    for snp in readVariantPositions(_in):
        snps.setdefault(snp[0], list()).append(VariantPosition(snp[0], int(snp[1]), snp[3], snp[4]))
    return snps

def processSNPs(_in, args):
    region_length = 2 * args.flanksize + 1

    with open('{}.filtered.tsv'.format(args.prefix), 'w') as filtered_out, \
         open('{}.unfiltered.tsv'.format(args.prefix), 'w') as unfiltered_out:

         for _id,_seq in readFasta(_in):
             for snp in snps.get(_id[1:], list()):
                 """
                 pull out flanking regions for each SNP
                 VCF is 1-based, so transform to base0
                 |left flank| = |right flank|, e.g. flanksize = 100 -> 2x100 (flanks) + 1 (SNP)
                 """
                 pos0 = snp.pos - 1
                 seq = _seq[pos0 - args.flanksize:pos0], _seq[pos0], _seq[pos0 + 1:pos0 + args.flanksize + 1]
                 assert snp.ref == seq[1]
                 """ blast """
                 query = '>{}:{}:{}{}\n{}{}{}\n'.format(snp.seqid, snp.pos, snp.ref, snp.alt, seq[0], seq[1], seq[2])
                 pr = sub.Popen(BLAST_CMD.format(args.blastdb), shell=True, stdin=sub.PIPE, stderr=sub.PIPE, stdout=sub.PIPE)
                 # print(BLAST_CMD.format(args.blastdb)) 
                 out, err = pr.communicate(query.encode()) 
                 out = out.decode()
                 print(out, file=unfiltered_out, flush=True)

                 if out.strip():
                     # print(out.strip().split('\n')[0].split('\t'), len(out.strip().split('\n')[0].split('\t')))
                     blast_hit = BlastHSP(*out.strip().split('\n')[0].split('\t'))
                     if filterBlastHit(blast_hit, region_length=region_length):
                         """
                         only filtered blast hits (qcov>=60%, ppos>=50%) are used for prediction,
                         the former criterion assures that the hit is spread over the whole region
                         and covers the SNP
                         """
                         predictEffect(seq, blast_hit, snp.ref, snp.alt, _out=filtered_out)



def main():
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

if __name__ == '__main__':
    main()
