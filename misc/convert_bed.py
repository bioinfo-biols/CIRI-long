#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from CIRI_long.align import GTFParser

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file, 'r') as f, open(out_file, 'w') as out:
    for line in f:
        if line.startswith('#'):
            continue
        content = line.rstrip().split('\t')
        parser = GTFParser(content)
        tmp_line = [parser.contig, parser.start, parser.end, parser.attr['circ_id'], 1_000, parser.strand, parser.start, parser.end]
        if parser.strand == "-":
            itemRgb = "43,140,190"
        else:
            itemRgb = "240,59,32"
        tmp_line.append(itemRgb)

        for iso in parser.attr['isoform'].split('|'):
            exons = iso.split(',')
            blockCount = len(exons)
            blockSize = []
            blockStarts = []
            for exon in exons:
                exon_st, exon_en = exon.split('-')
                blockSize.append(str(int(exon_en) -int(exon_st)))
                blockStarts.append(str(int(exon_st)-parser.start))
            out.write('\t'.join([str(x) for x in tmp_line + [blockCount, ','.join(blockSize), ','.join(blockStarts)]]) + '\n')