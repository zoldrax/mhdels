#!/usr/bin/env python3
import argparse
import sys
from collections import Counter
from math import log

import pysam


def entropy(s):
    p, lns = Counter(s), float(len(s))
    return sum(count / lns * -log(count / lns, 4) for count in p.values())


def comparefromstart(a, b):
    i = 0
    for x, y in zip(a.upper(), b.upper()):
        if x != y: break
        i += 1
    return (i)


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help="input VCF file", default=sys.stdin)
parser.add_argument('-o', '--output', help="output VCF file", default=sys.stdout)
parser.add_argument('-R', '--reference', help="refrence FASTA file", required=True)

args = parser.parse_args()

ref = pysam.FastaFile(args.reference)
vcf_in = pysam.VariantFile(args.input)
vcf_in.header.info.add('MHD_SEQ', '.', 'String', 'Microhomology sequence')
vcf_in.header.info.add('MHD_ENTLEN', 'A', 'Float', 'Entropy multiplyed by the length of microhomology sequence')
print(list(vcf_in.header.info))
vcf_out = pysam.VariantFile(args.output, 'w', header=vcf_in.header)
for record in vcf_in:
    start = record.stop
    stop = start + record.rlen - 1
    contig_len = ref.get_reference_length(record.contig)
    if stop > contig_len: stop = contig_len
    postfix = ref.fetch(reference=record.contig, start=start, end=stop)
    deletion = record.ref[len(record.alts[0]):]
    equal = comparefromstart(deletion, postfix)
    if (equal > 2):
        print(deletion)
        print(postfix)
        equalstr = deletion[:equal]
        print(equal * entropy(equalstr), equalstr)
        record.info['MHD_SEQ'] = equalstr
        record.info['MHD_ENTLEN'] = equal * entropy(equalstr)
    vcf_out.write(record)
vcf_out.close()
vcf_in.close()
ref.close()
