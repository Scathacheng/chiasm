#!/usr/bin/env python3

import gff
from argparse import ArgumentParser as AP
import logging as log
import random
from collections import defaultdict


cli = AP(description='''Generate mock circRNAs based on genome annotation''')
cli.add_argument('-g', '--gff', help='genome GFF annotation', metavar='FILE', required=True)
cli.add_argument('-n', '--number', type=int, default=100, metavar='N',
                 help='number of junctions (default: 100)')
cli.add_argument('-i', '--internal', action='store_true', help='exclude first and last exons')
cli.add_argument('-s', '--similar_length', metavar='FILE',
                 help='TAB separated file with 5\'-/3\'-flanking intron length')
cli.add_argument('-l', '--length', type=int, metavar='N', help='minimal length of flanking introns')
cli.add_argument('-v', '--verbose', action='store_true', help='print verbose output to STDERR')
args = cli.parse_args()
log.basicConfig(level=log.DEBUG if args.verbose else log.WARN, format='[%(levelname)s] %(message)s')

blur = 100


def weighted_choice(choices, weights):
    total = sum(weights)
    treshold = random.uniform(0, total)
    for k, weight in enumerate(weights):
        total -= weight
        if total < treshold:
            return choices[k]


def print_mock(exons, first, last, transcript):
    start = exons[first].start
    end = exons[last].end
    annot5 = genome.get_annotation_from_transcript(transcript, start)
    annot3 = genome.get_annotation_from_transcript(transcript, end)
    intron5 = genome.get_left_intron(exons[first])
    intron3 = genome.get_right_intron(exons[last])
    intron5 = 0 if intron5 is None else intron5.end-intron5.start
    intron3 = 0 if intron3 is None else intron3.end-intron3.start
    if not transcript.sense:
        intron5, intron3 = intron3, intron5
        annot5, annot3 = annot3, annot5
        first, last = len(exons)-last-1, len(exons)-first-1

    print('%s\t%d\t%d\tsplits:0:0:0:C:P\t0\t%c\t%s\t%d\t%d\t%s\t%s\t%d\t%d' %
          (transcript.chro, start, end, transcript.strand, transcript.name,
           first+1, last+1, annot5, annot3, intron5, intron3))


log.info('parsing genome annotation ...')
genome = gff.GFF(args.gff)

regions = {r.chro: r.end-r.start for _, r in genome.features.items() if r.type == 'region'}


log.info('generating %d random junctions ...' % args.number)
# print first line header
print('#ref\tstart\tend\tname\tscore\tstrand\ttranscript\t5exon\t3exon\t5annot\t3annot'
      '\t5intron\t3intron')


# only emit junctions with adjacent introns of comparable length to those in the file
if args.similar_length:
    introns = defaultdict(lambda: defaultdict(list))
    for t in [gene for name, gene in genome.features.items() if gene.type == 'mRNA']:
        exons = genome.get_exons_from_transcript(t, False)
        if len(exons) < 3:
            continue
        for i in range(1, len(exons)-1):
            for j in range(i, len(exons)-1):
                left = genome.get_left_intron(exons[i])
                left = int((left.end-left.start)/blur)
                right = genome.get_right_intron(exons[j])
                right = int((right.end-right.start)/blur)
                introns[left][right].append((exons, i, j, t))
    pairs = []
    with open(args.similar_length, 'r') as lines:
        for l in lines:
            first, second = [int(int(x)/blur) for x in l.rstrip().split('\t')[0:2]]
            pairs.append((first, second))
    log.info('sampling %d similar intron lengths' % len(pairs))

    if args.number < 1:
        args.number = len(pairs)

    count = 0
    while count <= args.number:
        s, e = random.choice(pairs)  # sample the number of junctions from the given distribution
        if s not in introns:
            log.warn('there is no upstream intron with approx. length %d' % (s*blur))
            continue
        if e not in introns[s]:
            log.warn('there is no right intron with approx. length %d' % (e*blur))
            continue
        if not introns[s][e]:
            log.info('too few candidates with intron lengths %d-%d' % (s*blur, e*blur))
            continue
        i = random.choice(range(len(introns[s][e])))
        exons, first, last, transcript = introns[s][e].pop(i)
        print_mock(exons, first, last, transcript)
        count += 1
    exit()


uniques = {}
first_index = 1 if args.internal else 0
if args.number < 1:
    for c in genome.genes:
        for s in genome.genes[c]:
            for g in genome.genes[c][s]:
                gene = genome.features[g]
                for t in genome.get_transcripts_from_gene(gene):
                    exons = genome.get_exons_from_transcript(t, False)
                    if args.internal and len(exons) < 3:
                        continue
                    last_index = len(exons)-1 if args.internal else len(exons)
                    for i in range(first_index, last_index):
                        for j in range(i, last_index):
                            loc = '%s|%d%c%d' % (t.chro, exons[i].start, t.strand, exons[j].end)
                            if loc in uniques:
                                continue
                            uniques[loc] = 1
                            if args.length and i > 0 and j < len(exons)-1:
                                if exons[i].start-exons[i-1].end < args.length or \
                                        exons[j+1].start-exons[j].end < args.length:
                                    continue
                            print_mock(exons, i, j, t)
    exit()


count = 0
while count < args.number:
    chro = weighted_choice(list(regions.keys()), weights=list(regions.values()))
    strand = random.choice(list(genome.genes[chro].keys()))
    genes = genome.genes[chro][strand]
    if not genes:
        continue
    gene = genome.features[random.choice(genes)]
    transcripts = genome.get_transcripts_from_gene(gene)
    if not transcripts:
        continue
    transcript = random.choice(transcripts)
    exons = genome.get_exons_from_transcript(transcript, False)    # sorted by pos
    if not exons:
        continue
    if args.internal and len(exons) < 3:
        continue
    last_index = len(exons)-2 if args.internal else len(exons)-1
    first = random.randint(first_index, last_index)
    last = random.randint(first, last_index)
    loc = '%s|%d%c%d' % (transcript.chro, exons[first].start, transcript.strand, exons[last].end)
    if loc in uniques:
        continue
    uniques[loc] = 1
    if args.length and first > 0 and last < len(exons)-1:
        if exons[first].start-exons[first-1].end < args.length or \
                exons[last+1].start-exons[last].end < args.length:
            continue
    log.info('found %d at %s' % (count, loc))
    print_mock(exons, first, last, transcript)
    count += 1
