#!/usr/bin/env python3
import sys
import gff
from argparse import ArgumentParser as AP
from Bio import SeqIO
import logging as log


cli = AP(description='Concatenates junctions based on BED format input and genome annotation')
cli.add_argument('-v', '--verbose', action='store_true', help='print progress info to STDERR')
cli.add_argument('-f', '--fasta', required=True, help='reference genome sequence', metavar='FILE')
cli.add_argument('-l', '--length', type=int, default=-1, metavar='N',
                 help=('output length of [exonic] sequence inside junction '
                       '(default: -1 for full seq)'))
cli.add_argument('-a', '--amend', type=int, default=0, metavar='N',
                 help='merge similar junctions [and correct by exon boundary] within N nt')
annot = cli.add_argument_group('annotation options')
annot.add_argument('-g', '--gff', help='GFF formated reference annotation', metavar='FILE')
annot.add_argument('-s', '--ignore_strand', action='store_true', help='ignore junction strand info')
annot.add_argument('-u', '--upstream', type=int, default=0, metavar='N',
                   help=('output 5`-flanking [intronic] sequence of length N'
                         '[-1 for complete intron] (default: 0)'))
annot.add_argument('-d', '--downstream', type=int, default=0, metavar='N',
                   help=('output 3`-flanking [intronic] sequence of length N '
                         '[-1 for complete intron] (default: 0)'))
cli.add_argument('--print_length_only', action='store_true',
                 help='output only flanking intron lengths')
args = cli.parse_args()
log.basicConfig(level=log.DEBUG if args.verbose else log.WARN, format='[%(levelname)s] %(message)s')
if args.length < 0 and not args.gff:
    log.warn('--length parameter can only be set to -1 if --gff annotation is specified. set to 0')
    args.length = 0


log.info('reading genome annotation form GFF/GTF...')
if args.gff:
    annot = gff.GFF(args.gff)


log.info('reading genome sequence from FASTA...')
refs = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))


# soft-limit up/downstream introns if values are negative
limit_upstream, limit_downstream = False, False
if args.upstream < -1:
    limit_upstream = True
    args.upstream = -args.upstream
if args.downstream < -1:
    limit_downstream = True
    args.downstream = -args.downstream

log.info('processing junctions from BED...')
count, exonic_count = 0, 0
for l in sys.stdin:
    l = l.rstrip()
    c = l.split('\t')
    if l.startswith('#') or len(c) < 6:
        continue
    count += 1
    chro, start, end, strand = c[0], int(c[1]), int(c[2]), c[5]
    exon_names = ''
    upstream, downstream = '', ''
    upstream = str(refs[chro].seq[start-1-args.upstream:start-1])
    downstream = str(refs[chro].seq[end:end+args.downstream])
    seq1 = refs[chro].seq[start-1:start-1+args.length]
    seq2 = refs[chro].seq[end-args.length:end]
    if strand == '-':
        downstream = str(refs[chro].seq[start-1-args.downstream:start-1].reverse_complement())
        upstream = str(refs[chro].seq[end:end+args.upstream].reverse_complement())
        seq1 = seq2.reverse_complement()
        seq2 = seq1.reverse_complement()

    if args.gff:
        if chro not in annot.genes:
            raise Exception('Reference "%s" is not part of the annotation' % chro)
            continue

        exons, trans = annot.get_exons(chro, start, end, None if args.ignore_strand else strand)

        if not exons:
            log.warn('%s|%d%c%d has no annotation. Skipping!' % (chro, start, strand, end))
            continue
        exonic_count += 1

        strand = exons[0].strand
        sense = strand == '+'
        exon_names = '%s %s' % (exons[0].name, exons[-1].name)
        if not sense:
            exon_names = '%s %s' % (exons[-1].name, exons[0].name)
        left = annot.get_left_intron(exons[0])
        right = annot.get_right_intron(exons[-1])
        if args.upstream == -1 or limit_upstream:
            if sense and not left or not sense and not right:
                log.warn('%s|%d%c%d has no upstream intron' % (chro, start, strand, end))
                upstream = ''
                continue
            if sense:
                s, e = left.start, left.end
                if limit_upstream and args.upstream < e-s:
                    s = e-args.upstream
                upstream = str(refs[chro].seq[s-1:e])
            else:
                s, e = right.start, right.end
                if limit_upstream and args.upstream < e-s:
                    e = s+args.upstream
                upstream = str(refs[chro].seq[s-1:e].reverse_complement())
        if args.downstream == -1 or limit_downstream:
            if sense and not right or not sense and not left:
                log.warn('%s|%d%c%d has no downstream intron' % (chro, start, strand, end))
                downstream = ''
                continue
            if sense:
                s, e = right.start, right.end
                if limit_downstream and args.downstream < e-s:
                    e = s+args.downstream
                downstream = str(refs[chro].seq[s-1:e])
            else:
                s, e = left.start, left.end
                if limit_downstream and args.downstream < e-s:
                    s = e-args.downstream
                downstream = str(refs[chro].seq[s-1:e].reverse_complement())

        # complete exonic length
        if args.length == -1:
            exon_names = ' '.join(sorted([x.name for x in exons]))
            seqs = []
            for x in exons:
                s = x.start if start < x.start else start   # print only sequence covered by junct
                e = x.end if x.end < end else end
                seqs.append(refs[chro].seq[s-1:e])
            seq1 = ''.join([str(s) for s in seqs] if sense else
                           [str(s.reverse_complement()) for s in reversed(seqs)])

    if args.print_length_only:
        print('%d\t%d' % (len(upstream), len(downstream)))
        continue

    print('>%s|%d%c%d %s\n%s%s%s%s'
          % (chro, start, strand, end, exon_names,
             upstream, str(seq1), str(seq2), downstream))

log.info('processed %d junctions (%d with exonic sequence)' % (count, exonic_count))
