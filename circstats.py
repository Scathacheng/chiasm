#!/usr/bin/env python3
import gff
# import re
import sys
import logging as log
from argparse import ArgumentParser as AP
# from collections import defaultdict

cli = AP(description='''Statistics of circRNAs based on genome annotation''')
cli.add_argument('-g', '--gff', help='genome GFF annotation', metavar='FILE')
cli.add_argument('-c', '--correct_strand', action='store_true',
                 help='correct strand information according to suitable transcript information')
cli.add_argument('-v', '--verbose', action='store_true',
                 help="print verbose output to STDERR")
args = cli.parse_args()
log.basicConfig(level=log.DEBUG if args.verbose else log.WARN, format='[%(levelname)s] %(message)s')


log.info('parsing annotation ...')
genome = gff.GFF(args.gff)


log.info('processing junctions ...')
num_junct, num_start, num_end, no_annot, num_corrected = 0, 0, 0, 0, 0
for l in sys.stdin:
    l = l.rstrip('\r\n')
    c = l.split('\t')
    if l.startswith('#') or len(c) < 6:
        print(l + '\tgene\ttranscript\t5exon\t3exon\t5annot\t3annot\t5intron\t3intron')
        continue
    chro, start, end, strand = c[0], int(c[1]), int(c[2]), c[5]

    if chro not in genome.genes:
        continue
    annot5, annot3 = '*', '*'
    parent_transcript, parent_gene = '*', '*'
    first_intron, last_intron = None, None
    exon5, exon3 = 0, 0
    # This is new
    gene, transcript, annot_s, annot_e = '', '', '*', '*'
    start_s, start_e, intron_s, intron_e = '', '', 0, 0
    num_s, num_e, best_s, best_e = 0, 0, '', ''
    transcripts = genome.get_transcripts(chro, start, end, None if args.correct_strand else strand)
    for t in transcripts:
        annot_s, annot_e = 'ncRNA', 'ncRNA'
        exons = genome.get_exons_from_transcript(t)
        transcript = t
        gene = t.parents if t.parents != t.name else gene
        for i, e in enumerate(exons):
            if e.start == start:
                best_s = e
                transcript = t
                num_s = i
            if e.end == end:
                best_e = e
                transcript = t
                num_e = i
            if best_s and best_e and best_s.parents == best_e.parents:
                if args.correct_strand:
                    strand = t.strand
                next_e, prev_e = num_e+1, num_s-1
                if strand == '+':
                    next_e, prev_e = prev_e, next_e
                intron_s = start-exons[prev_e].end if 0 < num_s else 0
                intron_e = exons[next_e].start-end if num_e < len(exons)-1 else 0
                break
        if best_s and best_e and best_s.parents == best_e.parents:
            break
    if transcript:
        # annot5 = '5UTR' if start <= cds.start else '3UTR' if cds.end < end else 'CDS'
        # annot3 = '5UTR' if start <= cds.start else '3UTR' if cds.end < end else 'CDS'
        annot_s = genome.get_annotation_from_transcript(transcript, start)
        annot_e = genome.get_annotation_from_transcript(transcript, end)
        gene = ','.join(transcript.parents)
        trans = transcript.name
    if strand == '-':
        intron_s, intron_e = intron_e, intron_s
        annot_s, annot_e = annot_e, annot_s
        num_s, num_e = num_e, num_s

    exons, transcript = genome.get_exons(chro, start, end, None if args.correct_strand else strand)
    if exons:
        if start == exons[0].start:
            num_start += 1
        if end == exons[-1].end:
            num_end += 1
        if start != exons[0].start and end != exons[-1].end:
            log.info('junction "%s|%d%c%d" is not exon boundary' % (chro, start, strand, end))
        parent_transcript = transcript.name
        parent_gene = transcript.parents[0]
        if args.correct_strand and strand != transcript.strand:
            num_corrected += 1
            strand = transcript.strand
        annot5 = genome.get_annotation_from_transcript(transcript, start)
        annot3 = genome.get_annotation_from_transcript(transcript, end)
        if annot5 == 'exonic':
            log.info('transcript %s has no CDS' % transcript.name)
        first_intron = genome.get_left_intron(exons[0])
        last_intron = genome.get_right_intron(exons[-1])
        all_exons = genome.get_exons_from_transcript(transcript, False)
        for i, e in enumerate(all_exons):
            if e.name == exons[0].name:
                exon5 = i+1  # if transcript.sense else len(all_exons)-i+1
            if e.name == exons[-1].name:
                exon3 = i+1  # if transcript.sense else len(all_exons)-i+1
    else:
        no_annot += 1
        log.info('junction "%s|%d%c%d" has no annotation' % (chro, start, strand, end))

    intron5 = first_intron.end-first_intron.start if first_intron else 0
    intron3 = last_intron.end-last_intron.start if last_intron else 0
    c[5] = strand
    if strand == '-':
        annot5, annot3 = annot3, annot5
        intron5, intron3 = intron3, intron5
        # exon5, exon3 = exon3, exon5

    num_junct += 1
    print('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d' %
          ('\t'.join(c), parent_gene, parent_transcript,
           exon5, exon3, annot5, annot3, intron5, intron3,
           gene, trans, num_s, num_e, annot_s, annot_e,
           intron_s, intron_e))

log.info('processed %d junctions (%d exon starts / %d exon ends / %d no annotation)' %
         (num_junct, num_start, num_end, no_annot))
if args.correct_strand:
    log.info('corrected strand information for %d junctions' % num_corrected)
