#!/usr/bin/env python3

from argparse import ArgumentParser as AP
import sys
import logging as log


cli = AP(description='Parses ODB relations and annotates a list with orthologs from a second')
cli.add_argument('-v', '--verbose', action='store_true', help='print verbose logging info')
cli.add_argument('-i', '--input_column', type=int, default=1, help='ID field/colum in input')
cli.add_argument('-o', '--ortho_column', type=int, default=1, help='ID field/colum in orthologs')
cli.add_argument('orthologs', help='list of orthologs with ID in first column')
args = cli.parse_args()
log.basicConfig(level=log.DEBUG if args.verbose else log.WARN, format='[%(levelname)s] %(message)s')


log.info('parsing orthologs ...')
ortho = {}
with open(args.orthologs) as f:
    for l in f:
        if l.startswith('#'):
            continue
        c = l.rstrip().split('\t')
        id = c[args.ortho_column-1]
        annot = '\t'.join([c[i] for i in range(len(c)) if i != args.ortho_column-1])
        if id in ortho:
            ortho[id].append(annot)
        else:
            ortho[id] = [annot]
log.info('%d orthologs ...' % len(ortho))


count, found = 0, 0
for l in sys.stdin:
    l = l.rstrip('\n\r ')
    if l.startswith('#'):
        print('%s\t%s' % (l, args.orthologs))
        continue
    count += 1
    id = l.split('\t')[args.input_column-1]
    add = ''
    if id in ortho:
        found += 1
        add = ','.join(ortho[id])
    print('%s\t%s' % (l, add))

log.info('found orthologs for %d in %d' % (found, count))
