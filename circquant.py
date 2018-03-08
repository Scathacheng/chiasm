#!/usr/bin/env python3

import re
import logging as log
# import gff
from collections import defaultdict
from argparse import ArgumentParser as AP


cli = AP(description='''Quantification of circRNA-Seq mapping results
Reads splice information from segemehls testrealign or circumjunct to quantify
circular junctions against their linear splicing counterparts and outputs a table
of all circular junctions found across multiple samples.''')
cli.add_argument('bedfiles', nargs='+', help='multiple BED files '
                 'from containing splice sites')
cli.add_argument('-c', '--min_circs', type=int, default=1, metavar='N',
                 help='minimum number of found circular junctions in total')
cli.add_argument('-e', '--min_experiments', type=int, default=1, metavar='N',
                 help='minimum number of samples containing a circular junction')
cli.add_argument('-i', '--min_insert', type=int, default=100, metavar='N',
                 help='minimal insert size between junctions (default=100)')
cli.add_argument('-I', '--max_insert', type=int, default=10000, metavar='N',
                 help='maximal insert size between junctions (default=10000)')
cli.add_argument('-a', '--amend', type=int, default=0, metavar='N',
                 help='merge similar junctions and correct by exon boundary within N nt')
cli.add_argument('-g', '--gff', help='genome GFF annotation', metavar='FILE')
cli.add_argument('-s', '--ignore_strand', action='store_true', help='ignore junction strand')
cli.add_argument('-S', '--correct_strand', action='store_true',
                 help='determine strand info from annotation')
norm = cli.add_mutually_exclusive_group()
norm.add_argument('-n', '--normalize', action='store_true',
                  help='normalize by total library junction reads')
norm.add_argument('-f', '--factors', type=float, nargs='+', metavar='F',
                  help='normalize using these factors')
norm.add_argument('-H', '--host_counts', metavar='FILE',
                  help='normalize by host gene featureCounts in each library')
cli.add_argument('-v', '--verbose', action='store_true',
                 help="print verbose output to STDERR")
args = cli.parse_args()
log.basicConfig(level=log.DEBUG if args.verbose else log.WARN, format='[%(levelname)s] %(message)s')


junctions, splice, counts = defaultdict(int), {}, {}  # global junctions, L/R linear, circ per file
if args.factors and len(args.factors) != len(args.bedfiles):
    log.error('number of normalization factors does not match number of BED files: %d != %d' %
              (len(args.factors), len(args.bedfiles)))
    exit(1)

if args.host_counts:
    exp = '\t'.join(['%s_junct\thost_FPKM\tnormalized' % bed for bed in args.bedfiles])
    with open(args.host_counts, 'r') as countfile:
        header = countfile.readline()
        header = countfile.readline().rstrip()
        c = header.split('\t')
        if len(c)-6 != len(args.bedfiles):
            log.warn('number of junction files deviates from number of BAM files in %s: %d != %d' %
                     (args.host_counts, len(args.bedfiles), len(c)-6))
            exit(1)

        for i, l in enumerate(c[6:]):
            log.info('matching BED: %s\t to BAM: %s' % (args.bedfiles[i], l))

        libs = len(c)-6
        for line in countfile:
            c = line.rstrip().split('\t')
            gene, length = c[0], int(c[5])
            counts[gene] = [int(c[6+i])*1000/length if length > 0 else 0 for i in range(libs)]


log.info('reading junctions from %d files ...' % len(args.bedfiles))
for bed in args.bedfiles:     # this script can be called with multiple BED files
    splice[bed] = {'total': 0}
    with open(bed, 'r') as bedfile:
        for line in bedfile:
            c = line.rstrip().split()
            if line.startswith('#') or len(c) < 6:
                continue
            chro, start, end, splits, strand = c[0], int(c[1]), int(c[2]), c[3], c[5]
            if args.ignore_strand:
                strand = '+'
            s = splits.split(':')               # splits:12:12:14:C:P
            typ, flag, passed = s[0], s[4], s[5]
            junct, left, right = int(s[1]), int(s[2]), int(s[3])
            splice[bed]['total'] += junct

            # circular split which passed or failed
            if typ == 'splits' and flag == 'C' and passed in ['P', 'F']:
                key = '%s\t%d\t%c\t%d' % (chro, start, strand, end)
                splice[bed][key] = [left, junct, right]
                junctions[key] += junct
log.info('parsed %d individual junctions' % len(junctions))


# merge junctions
if args.amend:
    merge_count = 0
    for j in sorted(junctions.keys()):    # iterate over a copy of the keys
        merged = False
        for a in range(-args.amend, args.amend+1):      # from -x to x is easier than 0 to x for
            for b in range(-args.amend, args.amend+1):  # positive and negative, take any, not best
                if a == 0 and b == 0:
                    continue
                chro, start, strand, end = j.split('\t')
                start, end = int(start), int(end)
                p = '%s\t%d\t%c\t%d' % (chro, start+a, strand, end+b)
                if p in junctions:
                    if junctions[j] < junctions[p]:     # merge weaker into stronger junctions
                        merged = True
                        merge_count += 1
                        junctions[p] += junctions.pop(j)
                        for bed in splice:
                            if j in splice[bed]:
                                s = splice[bed].pop(j)
                                if p not in splice[bed]:
                                    splice[bed][p] = [0, 0, 0]
                                splice[bed][p][0] += s[0]
                                splice[bed][p][1] += s[1]
                                splice[bed][p][2] += s[2]
                        break
            if merged:
                break
    log.info('merged %d junctions ...' % merge_count)


overwritten = 0
starts, ends, exons, features, strands = {}, {}, {}, {}, {}
if args.gff:
    log.info('reading GFF annotation file ...')
    # genome = gff.GFF(args.gff)
    with open(args.gff, 'r') as gff:
        for line in gff:
            c = line.split('\t')
            if line.startswith('#') or len(c) < 9:
                continue
            chro, typ, start, end, strand = c[0], c[2], int(c[3]), int(c[4]), c[6]

            m = re.search('ID=([^;]+);?', line)
            if not m:
                # log.warn('did not find ID')
                continue
            name = m.group(1)

            gene = ''
            m = re.search('gene=([^;]+);?', line)
            if m:
                gene = m.group(1)

            parent = ''
            m = re.search('Parent=([^;]+);?', line)
            if m:
                parent = m.group(1)

            features[name] = {'chro': chro, 'typ': typ, 'start': start,
                              'end': end, 'strand': strand, 'gene': gene, 'parent': parent}

            strands[name] = strand
            if args.ignore_strand:
                strand = '+'

            if typ == 'exon':
                if chro not in exons:  # initialize dicts
                    starts[chro] = {'+': {}, '-': {}}
                    ends[chro] = {'+': {}, '-': {}}
                    exons[chro] = {'+': {}, '-': {}}

                if start in exons[chro][strand]:    # this is dirty and can lead to problems
                    overwritten += 1
                starts[chro][strand][start] = name
                ends[chro][strand][end] = name
                exons[chro][strand][start] = (end, name)
    log.info('read %d exons (%d redundant) ...' %
             (len([e for c in exons for s in exons[c] for e in exons[c][s]]), overwritten))


exp = '\t'.join(['%s_junct\tleft\tright' % bed for bed in args.bedfiles])
if args.factors:
    exp = '\t'.join(['%s_junct\tnormalized' % bed for bed in args.bedfiles])
    for i, f in enumerate(args.bedfiles):
        splice[f]['factor'] = args.factors[i]
    args.normalize = True
elif args.normalize:
    log.info('using normalization factors: library / junctions / factor')
    exp = '\t'.join(['%s_junct\tnormalized' % bed for bed in args.bedfiles])
    avg = sum([splice[f]['total'] for f in splice])/len(splice)
    for f in splice:
        splice[f]['factor'] = avg/splice[f]['total']
        log.info('%s\t%d\t%.4f' % (f, splice[f]['total'], splice[f]['factor']))
elif args.host_counts:
    exp = '\t'.join(['%s_junct\thost_FPKM\tnormalized' % bed for bed in args.bedfiles])


# print first line header
print('#ref\tstart\tend\tname\tscore\tstrand\tsum_linear\tleft_exon\tright_exon\t'
      'parent_gene\t%s' % exp)

log.info('evaluating insert length ...')
start_found, end_found, inexact_start, inexact_end = 0, 0, 0, 0
insert_fail, sample_fail, circ_fail, num_corrected = 0, 0, 0, 0
for j in sorted(junctions.keys()):
    if junctions[j] < args.min_circs:
        circ_fail += 1
        continue
    chro, start, strand, end = j.split('\t')
    start, end = int(start), int(end)

    # amend to annotation and add info
    start_exon, end_exon, insert = '*', '*', end-start+1
    start_parent, end_parent, gene = '*', '*', '*'
    num_exons, insert, before, after, introns, = 0, 0, '', '', []
    if chro in exons:
        # if chro in genome.genes:

        # align with nearby exon boundaries
        sx = starts[chro][strand]
        ex = ends[chro][strand]
        for i in range(-args.amend, args.amend+1):      # from -x to x is easier than 0 to x for
            if start_exon == '*' and start+i in sx:
                if i != 0:
                    inexact_start += 1
                    start += i
                start_exon = sx[start]
                start_found += 1
            if end_exon == '*' and end+i in ex:
                if i != 0:
                    inexact_end += 1
                    end += i
                end_exon = ex[end]
                end_found += 1

        # get parent gene
        if start_exon in features:
            start_parent = features[start_exon]['parent']
        if end_exon in features:
            end_parent = features[end_exon]['parent']

        gene = 'intergenic'
        if start_parent != '*' and end_parent != '*':
            for sp in start_parent.split(','):
                for ep in end_parent.split(','):
                    if sp == ep and sp in features:
                        gene = features[sp]['parent'] if features[sp]['parent'] else sp
                        continue
                    if sp in features and ep in features:
                        sg, eg = features[sp]['parent'], features[ep]['parent']
                        if sg == eg:
                            gene = sg
                            continue

        # compute exonic and intronic length
        last = start-1
        for s in sorted(exons[chro][strand]):
            e, n = exons[chro][strand][s]
            if e < start:   # exons before junction
                before = n
                continue
            if end < s:     # exon after junction
                after = n
                break
            if end < e:     # exons just ending after junction
                continue
            if last < s:
                num_exons += 1
                if start-1 < last:
                    introns.append(s-last)
                insert += e-s+1
                last = e

    # we really only care about a single gene or the info, that more than one gene is involved
    if args.host_counts and gene not in counts:     # with no count info, we cannot normalize
        if gene not in ['*', 'intergenic']:
            log.warn('Did not find "%s" in host counts!' % gene)
        continue

    if args.correct_strand and gene in strands:
        if strand != strands[gene]:
            num_corrected += 1
        strand = strands[gene]

    start_intron = start-features[before]['end'] if before else 0
    end_intron = features[after]['start']-end if after else 0
    intron_length = sum(introns)/len(introns) if introns else 0.0

    if not args.min_insert <= insert < args.max_insert:
        insert_fail += 1
        continue

    lefts = sum([splice[f][j][0] for f in splice if j in splice[f]])
    rights = sum([splice[f][j][2] for f in splice if j in splice[f]])
    exp = ''
    sample_circ, lin = 0, 0     # initialize counters
    for i, bed in enumerate(args.bedfiles):     # columns for each BED file
        left, junct, right = 0, 0, 0
        if j in splice[bed]:
            left, junct, right = splice[bed][j]
            sample_circ += 1    # circ junctions in this sample
        lin += left+right
        if args.normalize:
            factor = splice[bed]['factor']
            exp += '\t%d\t%.2f' % (junct, junct*factor)
        elif args.host_counts:
            fpkm = counts[gene][i]
            exp += '\t%d\t%.2f\t%.5f' % (junct, fpkm, junct/fpkm if fpkm > 0 else 0)
        else:
            exp += '\t%d\t%d\t%d' % (junct, left, right)
    text = ('%s\t%d\t%d\tcirc:%d:%d:%d:C:P\t%d\t%c\t%d\t%s\t%s\t%s%s' %
            (chro, start, end, junctions[j], lefts, rights, junctions[j], strand,
             lin, start_exon, end_exon, gene, exp))
    # print if min. total junctions and min. number of samples with circs is met
    if sample_circ < args.min_experiments:
        sample_fail += 1
        continue
    print(text)

log.info('removed %d junctions (insert: %d, few reads: %d, few samples: %d)' %
         (insert_fail+circ_fail+sample_fail, insert_fail, circ_fail, sample_fail))

if args.gff:
    log.info('exonic junction start: %d (%d amended) / exonic junction ends: %d (%d amended)' %
             (start_found, inexact_start, end_found, inexact_end))
    if args.correct_strand:
        log.info('corrected strand information for %d junctions' % num_corrected)
