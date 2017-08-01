#!/usr/bin/env python
# coding=utf-8
import sys
import csv
import subprocess as sub
import os

from ktoolu_io import readFasta
from extractFeatures import translate

BLAST_DB = os.path.join(os.path.dirname(sys.argv[0]), 'grass_proteins.cd95.fa')
BLAST_CMD = 'blastx -max_target_seqs 1 -db %s -outfmt "6 std sseqid qstart qend sstart send evalue bitscore qlen slen positive gaps ppos qframe staxids salltitles qseq sseq"' % BLAST_DB
MARKER_LENGTH = 201
"""
query_170686	EMT30199.1	81.818	66	12	0	199	2	385	450	1.76e-19	83.2	gi|475617667|gb|EMT30199.1|	201	1052	59	0	89.39	-3	N/A	Putative disease resistance protein RGA3 [Aegilops tauschii]	DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM	DDIFHEMLKDITGDRHSHISDHEELEEKLKKELHGKRFFLILDDLWVKTKNDPQLEELISPLNVGM

FILTERING ['query_170686', 'EMT30199.1', '81.818', '66', '12', '0', '199', '2', '385', '450', '1.76e-19', '83.2', 'gi|475617667|gb|EMT30199.1|', '201', '1052', '59', '0', '89.39', '-3', 'N/A', 'Putative disease resistance protein RGA3 [Aegilops tauschii]', 'DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM', 'DDIFHEMLKDITGDRHSHISDHEELEEKLKKELHGKRFFLILDDLWVKTKNDPQLEELISPLNVGM']
True True
['query_170686', 'EMT30199.1', '81.818', '66', '12', '0', '199', '2', '385', '450', '1.76e-19', '83.2', 'gi|475617667|gb|EMT30199.1|', '201', '1052', '59', '0', '89.39', '-3', 'N/A', 'Putative disease resistance protein RGA3 [Aegilops tauschii]', 'DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM', 'DDIFHEMLKDITGDRHSHISDHEELEEKLKKELHGKRFFLILDDLWVKTKNDPQLEELISPLNVGM'] TCATCCCAACACTGAGCGGAGAGATTACTTCTACCAGCTGTGGGTCATTCTTGTTCTTCACCCAGAGATCATCCAATATCAAGAAGAAACGTTTGCCACGCAATGATTTCTTCAGCTTCTCTTCCAGCTCCTCACGATCTGAAATATAGGAGTGCCGATTTTTGGTAATATCCTTCAGCATTTCATTAAATATAGCATCCA G A
FILTERING ['query_719213', 'EMT33476.1', '71.212', '66', '19', '0', '3', '200', '756', '821', '1.08e-26', '103', 'gi|475626793|gb|EMT33476.1|', '201', '1051', '53', '0', '80.30', '3', 'N/A', 'Disease resistance RPP8-like protein 3 [Aegilops tauschii]', 'FPLLSFMNISVDRVQPEVDIQILGMFPALRVLRLWANKRRYSCVEMFVVGANAFPCLRECLFSGFL', 'FPLLSSMSIEVDRVRPEVDIQILGKLPALRFLWLWVNKSQHTRVETFVIGANAFPCLRECRFHEFL']
True True
['query_719213', 'EMT33476.1', '71.212', '66', '19', '0', '3', '200', '756', '821', '1.08e-26', '103', 'gi|475626793|gb|EMT33476.1|', '201', '1051', '53', '0', '80.30', '3', 'N/A', 'Disease resistance RPP8-like protein 3 [Aegilops tauschii]', 'FPLLSFMNISVDRVQPEVDIQILGMFPALRVLRLWANKRRYSCVEMFVVGANAFPCLRECLFSGFL', 'FPLLSSMSIEVDRVRPEVDIQILGKLPALRFLWLWVNKSQHTRVETFVIGANAFPCLRECRFHEFL'] TGTTTCCCCTCCTCTCCTTCATGAATATTTCAGTGGACAGAGTCCAGCCTGAAGTCGACATTCAGATCCTCGGGATGTTTCCTGCTCTTCGTGTTCTCAGGCTTTGGGCCAACAAACGTAGGTACAGTTGTGTCGAAATGTTCGTTGTTGGCGCTAATGCATTCCCGTGTTTGAGAGAGTGCCTCTTCTCTGGTTTTCTCA G A
FILTERING ['query_719427', 'EMS54357.1', '77.273', '66', '15', '0', '2', '199', '671', '736', '6.85e-29', '110', 'gi|474071839|gb|EMS54357.1|', '201', '831', '56', '0', '84.85', '2', 'N/A', 'Disease resistance protein RPP13 [Triticum urartu]', 'FPRGAMPMLEILWFHARASDIAGGDLDVSMGHLPSLQQVQVGFWREEGSSSDKCKEDADVLLRHAA', 'FPRRAMSRLEILRFQARSSDIASGDLDVGMGHLPSLQEVRVGLWLEKGSSSDKCKGDADVVLRHAA']
True True
['query_719427', 'EMS54357.1', '77.273', '66', '15', '0', '2', '199', '671', '736', '6.85e-29', '110', 'gi|474071839|gb|EMS54357.1|', '201', '831', '56', '0', '84.85', '2', 'N/A', 'Disease resistance protein RPP13 [Triticum urartu]', 'FPRGAMPMLEILWFHARASDIAGGDLDVSMGHLPSLQQVQVGFWREEGSSSDKCKEDADVLLRHAA', 'FPRRAMSRLEILRFQARSSDIASGDLDVGMGHLPSLQEVRVGLWLEKGSSSDKCKGDADVVLRHAA'] GTTTCCACGAGGAGCTATGCCAATGCTTGAAATCCTGTGGTTCCATGCCCGGGCGTCGGATATTGCCGGCGGTGACCTGGATGTCAGCATGGGGCACCTCCCTTCCCTCCAGCAAGTCCAGGTTGGCTTTTGGCGTGAGGAAGGCAGCTCTTCAGATAAGTGCAAGGAGGACGCAGATGTCCTGCTGAGGCATGCAGCGGA C T
FILTERING ['query_723700', 'EMT21329.1', '98.507', '67', '1', '0', '201', '1', '464', '530', '7.59e-37', '132', 'gi|475589981|gb|EMT21329.1|', '201', '1467', '66', '0', '98.51', '-1', 'N/A', 'Putative disease resistance protein RGA4 [Aegilops tauschii]', 'LSSISDHRGLNKKLKEALRGKRFLLILDDLWVKNKNDQQLEELISPLNVGLKGSKILVTARTKEAAG', 'LPSISDHRGLNKKLKEALRGKRFLLILDDLWVKNKNDQQLEELISPLNVGLKGSKILVTARTKEAAG']

GGATGCTATATTTAATGAAATGCTGAAGGATATTACCAAAAATCGGCACTCCTATATTTCAGATCGTGAGGAGCTGGAAGAGAAGCTGAAGAAATCATTGCGTGGCAAACGTTTCTTCTTGATATTGGATGATCTCTGGGTGAAGAACAAGAATGACCCACAGCTGGTAGAAGTAATCTCTCCGCTCAGTGTTGGGATGAA
TCATCCCAACACTGAGCGGAGAGATTACTTCTACCAGCTGTGGGTCATTCTTGTTCTTCACCCAGAGATCATCCAATATCAAGAAGAAACGTTTGCCACG C AATGATTTCTTCAGCTTCTCTTCCAGCTCCTCACGATCTGAAATATAGGAGTGCCGATTTTTGGTAATATCCTTCAGCATTTCATTAAATATAGCATCCA
marker[1-201]: TCATCCCAACACTGAGCGGAGAGATTACTTCTACCAGCTGTGGGTCATTCTTGTTCTTCACCCAGAGATCATCCAATATCAAGAAGAAACGTTTGCCACGCAATGATTTCTTCAGCTTCTCTTCCAGCTCCTCACGATCTGAAATATAGGAGTGCCGATTTTTGGTAATATCCTTCAGCATTTCATTAAATATAGCATCCA

marker[2-199]: CATCCCAACACTGAGCGGAGAGATTACTTCTACCAGCTGTGGGTCATTCTTGTTCTTCACCCAGAGATCATCCAATATCAAGAAGAAACGTTTGCCACGCAATGATTTCTTCAGCTTCTCTTCCAGCTCCTCACGATCTGAAATATAGGAGTGCCGATTTTTGGTAATATCCTTCAGCATTTCATTAAATATAGCATC
marker_rc[2-199]: GATGCTATATTTAATGAAATGCTGAAGGATATTACCAAAAATCGGCACTCCTATATTTCAGATCGTGAGGAGCTGGAAGAGAAGCTGAAGAAATCATTGCGTGGCAAACGTTTCTTCTTGATATTGGATGATCTCTGGGTGAAGAACAAGAATGACCCACAGCTGGTAGAAGTAATCTCTCCGCTCAGTGTTGGGATG

blastxseq: DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM
marker_rc[2-199]_aa:           DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM (identical)

['query_170686', 'EMT30199.1', '81.818', '66', '12', '0', '200', '3', '385', '450', '1.76e-19', '83.2', 'gi|475617667|gb|EMT30199.1|', '201', '1052', '59', '0', '89.39', '-2', 'N/A', 'Putative disease resistance protein RGA3 [Aegilops tauschii]', 'DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM', 'DDIFHEMLKDITGDRHSHISDHEELEEKLKKELHGKRFFLILDDLWVKTKNDPQLEELISPLNVGM'] 
GGATGCTATATTTAATGAAATGCTGAAGGATATTACCAAAAATCGGCACTCCTATATTTCAGATCGTGAGGAGCTGGAAGAGAAGCTGAAGAAATCATTGCGTGGCAAACGTTTCTTCTTGATATTGGATGATCTCTGGGTGAAGAACAAGAATGACCCACAGCTGGTAGAAGTAATCTCTCCGCTCAGTGTTGGGATGAA
G T C 
GGATGCTATATTTAATGAAATGCTGAAGGATATTACCAAAAATCGGCACTCCTATATTTCAGATCGTGAGGAGCTGGAAGAGAAGCTGAAGAAATCATTGTTGGCAAACGTTTCTTCTTGATATTGGATGATCTCTGGGTGAAGAACAAGAATGACCCACAGCTGGTAGAAGTAATCTCTCCGCTCAGTGTTGGGATGAA GCYI**NAEGYYQKSALLYFRS*GAGREAEEIIVGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM -2 1


['query_170686', 'EMT30199.1', '81.818', '66', '12', '0', '200', '3', '385', '450', '1.76e-19', '83.2', 'gi|475617667|gb|EMT30199.1|', '201', '1052', '59', '0', '89.39', '-2', 'N/A', 'Putative disease resistance protein RGA3 [Aegilops tauschii]', 'DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM', 'DDIFHEMLKDITGDRHSHISDHEELEEKLKKELHGKRFFLILDDLWVKTKNDPQLEELISPLNVGM']
TTCATCCCAACACTGAGCGGAGAGATTACTTCTACCAGCTGTGGGTCATTCTTGTTCTTCACCCAGAGATCATCCAATATCAAGAAGAAACGTTTGCCAC-G-CAATGATTTCTTCAGCTTCTCTTCCAGCTCCTCACGATCTGAAATATAGGAGTGCCGATTTTTGGTAATATCCTTCAGCATTTCATTAAATATAGCATCC
G A 
GGATGCTATATTTAATGAAATGCTGAAGGATATTACCAAAAATCGGCACTCCTATATTTCAGATCGTGAGGAGCTGGAAGAGAAGCTGAAGAAATCATTG-T-GTGGCAAACGTTTCTTCTTGATATTGGATGATCTCTGGGTGAAGAACAAGAATGACCCACAGCTGGTAGAAGTAATCTCTCCGCTCAGTGTTGGGATGAA DAIFNEMLKDITKNRHSYISDREELEEKLKKSLCGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM -2 1

DAIFNEMLKDITKNRHSYISDREELEEKLKKSLRGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM
DAIFNEMLKDITKNRHSYISDREELEEKLKKSLCGKRFFLILDDLWVKNKNDPQLVEVISPLSVGM
"""

def reverseComplement(seq):
    """
    Returns the reverse complement of nucleic acid seqence input.
    """
    compl= dict(zip('ACGTNRYWSMKBHDV', 'TGCANYRWSKMVDHB'))
    return ''.join([compl[base]
                    for base in seq.upper().replace('U', 'T')])[::-1]


def filterBlastHit(blast_hit):
    # print 'FILTERING', blast_hit
    # print float(blast_hit[2]) >= 50, (3 * float(blast_hit[3])) / MARKER_LENGTH >= 0.75
    return float(blast_hit[17]) >= 50 and (3 * float(blast_hit[3])) / MARKER_LENGTH >= 0.5 # 0.75
   
def predictEffect(seq, blast_hit, ref_allele, alt_allele, _out): 
    def extractCodon(seq, strand='+'):
        cpos = len(seq[0]) % 3
        if cpos == 0:
            return seq[1] + seq[2][:2], cpos
        elif cpos == 1:
            return seq[0][-1] + seq[1] + seq[2][0], cpos
        elif cpos == 2:
            return seq[0][-2:] + seq[1], cpos
       
    frame = int(blast_hit[18])
    qstart, qend = sorted(map(int, blast_hit[6:8]))
    varseq = seq[0] + seq[1] + seq[2]
    if qend < qstart:
        qstart, qend = qend, qstart
        varseq = reverseComplement(seq[2]) + reverseComplement(seq[1]) + reverseComplement(seq[0])

    seq, strand = (seq[0][qstart-1:], seq[1], seq[2][:qend-len(seq[0])-1]), '+'
    if frame < 0:
        seq, strand = tuple(reversed(map(reverseComplement, seq))), '-'
        ref_allele, alt_allele = map(reverseComplement, [ref_allele, alt_allele])
 
    codon, cpos = extractCodon(seq, strand)
    alt_codon = list(codon)
    alt_codon[cpos] = alt_allele
    alt_codon = ''.join(alt_codon)
     
        
    varseq = seq[0] + seq[1] + seq[2]  
    # print blast_hit, seq, ref_allele, alt_allele, varseq, translate(varseq), frame, abs(frame)-1, codon, translate(codon), cpos, alt_codon, translate(alt_codon)
    out = [blast_hit[0], blast_hit[12], blast_hit[20], blast_hit[2], abs(float(blast_hit[7])-float(blast_hit[6]))/float(blast_hit[13]), blast_hit[10]]
    mclass = 'nonsense' if translate(alt_codon) == '*' else ('missense' if translate(codon) != translate(alt_codon) else 'silent')
    out.extend([codon, translate(codon), alt_codon, translate(alt_codon), mclass])
    _out.write('\t'.join(map(str, out)) + '\n')
    _out.flush()
    


def readVariantPositions(_in):
    for row in csv.reader(_in, delimiter='\t'):
        if not row[0].startswith('#'):
            yield row

def readSNPs(_in):
    snps = dict()
    for snp in readVariantPositions(sys.stdin): 
        if snp[0] not in snps:
            snps[snp[0]] = list()
        snps[snp[0]].append((int(snp[1]), snp[3], snp[4]))
    return snps

def processSNPs(_in):
   with open('scvep.filtered.tsv', 'w') as filtered_out, open('scvep.unfiltered.tsv', 'w') as unfiltered_out:
       for _id,_seq in readFasta(_in):
           for pos, ref, alt in snps.get(_id[1:], list()):
               pos0 = pos - 1
               seq = _seq[pos0 - 100:pos0], _seq[pos0], _seq[pos0 + 1:pos0 + 101]
               # unfiltered_out.write('\t'.join((_id, ''.join(seq))) + '\n')
               pr = sub.Popen(BLAST_CMD, shell=True, stdin=sub.PIPE, stderr=sub.PIPE, stdout=sub.PIPE)
               out, err = pr.communicate('>query_%s\n%s%s%s\n' % (pos, seq[0], seq[1], seq[2]))
               unfiltered_out.write(out) #, out.split('\t')[17]
               unfiltered_out.flush()
               if out.strip():
                   blast_hit = out.strip().split('\n')[0].split('\t')
                   if filterBlastHit(blast_hit):
                       # print('PASS')
                       predictEffect(seq, blast_hit, ref, alt, _out=filtered_out)

if __name__ == '__main__':
    snps = readSNPs(sys.stdin)        
    with open(sys.argv[1]) as seq_in:
        processSNPs(seq_in) 
         
