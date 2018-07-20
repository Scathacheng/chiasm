#!/usr/bin/env python3

import regex
from subprocess import Popen, PIPE
from Bio import SeqIO
import argparse
from sys import stderr as err
import logging as log

description = 'Complementary Oligonucleotide Interaction Finder For RNA\n' \
    'Searches for regions in the reference FASTA that are complementary to the query sequence. ' \
    'Allowing wobble base pairs between G:U, insertions, deletions and substitutions (errors), ' \
    'if set in the parameters'
cli = argparse.ArgumentParser(prog='coiffR', description=description)
cli.add_argument('-f', '--fasta', type=argparse.FileType(), required=True,
                 help='FASTA file of reference genome')
cli.add_argument('-q', '--query', required=True, metavar='SEQ', help='oligonucleotide sequence')
cli.add_argument('-r', '--revcomp', action='store_true', help='also search in reverse '
                 'complement of reference')
cli.add_argument('-e', '--errors', type=int, default=0, metavar='N', help='max. allowed errors')
cli.add_argument('-i', '--inserts', type=int, default=0, metavar='N', help='max. allowed inserts')
cli.add_argument('-d', '--deletions', type=int, default=0, metavar='N',
                 help='max. allowed deletions')
cli.add_argument('-s', '--substitutions', type=int, default=0, metavar='N',
                 help='max. allowed substitutions')
cli.add_argument('-w', '--wobble', action='store_true', help='allow wobble base pairs')
cli.add_argument('-m', '--mirna', type=int, default=0, metavar='N',
                 help='miRNA target prediction mode with min. N matches')
cli.add_argument('-F', '--fold', action='store_true', help='include free energy score')
cli.add_argument('-v', '--verbose', action='store_true', help='print verbose output to STDERR')
args = cli.parse_args()
log.basicConfig(level=log.DEBUG if args.verbose else log.WARN, format='[%(levelname)s] %(message)s')

TR = {'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A', 'N': 'N'}
WOBBLE = {'G': 'CT', 'T': 'AG'}
MIRNA_SEED = 6+1  # add 1 for first position usually not matching


# reverse complements the query before searching!
def query_pattern(query, errors=0, inserts=0, deletions=0, substitutions=0, wobble=False):
    tr = TR.copy()
    tr['N'] = '[ACGT]'
    if wobble:
        tr['G'] = '[CT]'
        tr['T'] = '[AG]'
    p = ''
    for c in query[::-1]:
        if c not in tr:
            p += c
            err.write('encountered illegal symbol in query: %s' % c)
        else:
            p += tr[c]
    fuzzy = []
    if inserts > 0:
        fuzzy.append('i<=%d' % inserts)
    if deletions > 0:
        fuzzy.append('d<=%d' % deletions)
    if substitutions > 0:
        fuzzy.append('s<=%d' % substitutions)
    if errors > 0:
        fuzzy.append('e<=%d' % errors)
    if fuzzy:
        p = '(%s){%s}' % (p, ','.join(fuzzy))
    log.info('searching with pattern: "%s"' % p)
    return regex.compile(p)


def mirna_pattern(query):
    tr = TR.copy()
    tr['N'] = '[ACGT]'
    if args.wobble:
        tr['G'] = '[CT]'
        tr['T'] = '[AG]'
    seed = MIRNA_SEED if MIRNA_SEED < len(query) else len(query)
    p = '.{%d}' % (len(query)-seed)
    for c in query[seed:1:-1]:
        if c not in tr:
            raise Exception('encountered illegal symbol in query: %s' % c)
        p += tr[c]
    p += '.'
    return regex.compile(p)


def score(x, y):
    if x == '-' or y == '-':
        return -1
    if y in TR[x]:
        return 1
    if x in WOBBLE and y in WOBBLE[x] or y in WOBBLE and x in WOBBLE[y]:
        return 0.5
    return -0.5


def align(a, b):
    errors = args.errors + args.inserts + args.deletions + args.substitutions
    if not args.mirna and not args.wobble and errors == 0:
        return (a, b, ''.join(['|' for _ in range(len(a))]), len(a))
    m = len(a)
    n = len(b)
    H = [[0]*(n+1) for _ in range(m+1)]
    T = [[0]*(n+1) for _ in range(m+1)]
    C = [[' ']*(n+1) for _ in range(m+1)]
    best = 0
    for i in range(m):  # fill matrices
        for j in range(n):
            match = H[i][j] + score(a[i], b[j])  # match
            delete = H[i+1][j] + score('-', b[j])     # deletion
            insert = H[i][j+1] + score(a[i], '-')     # insertion
            h = max(match, delete, insert)
            H[i+1][j+1] = h
            if h > best:
                best = h
            if h == match:
                T[i][j] = 1
                if a[i] == TR[b[j]]:
                    C[i][j] = '|'
                elif a[i] in WOBBLE and b[j] in WOBBLE[a[i]]:
                    C[i][j] = ':'
                else:
                    C[i][j] = '*'
            elif h == delete:
                T[i][j] = 2
            else:
                T[i][j] = 3
    i = m-1
    j = n-1
    ra = ''
    rb = ''
    cig = ''
    while i >= 0 and j >= 0:    # backtrack
        t = T[i][j]
        cig += C[i][j]
        if t == 1:
            ra += a[i]
            rb += b[j]
            i -= 1
            j -= 1
        elif t == 3:
            ra += a[i]
            rb += '-'
            i -= 1
        elif t == 2:
            ra += '-'
            rb += b[j]
            j -= 1
        else:
            break
    while i >= 0:   # fill b with '-'
        i -= 1
        ra += a[i]
        rb += '-'
        cig += ' '
    # while j >= 0:   # fill a with '-'
    #     j -= 1
    #     ra += '-'
    #     rb += b[j]
    #     cig += ' '
    return (ra, rb, cig, best)


def fold(a, b):
    r = Popen('echo "%s&%s" | RNAcofold --noPS' % (a, b), shell=True, stdout=PIPE).communicate()[0]
    if not r:
        return 0.0
    m = regex.search(r'\(\s*(\-?\d+\.\d+)\)', r.decode('utf-8'))
    if not m:
        return 0.0
    return float(m.group(1))


def mirna_score(cigar):
    bit = 0
    for i, c in enumerate(cigar[::-1]):
        if c == '|':
            bit += 1
        else:
            if 0 < i < 9:
                if c == ':':
                    bit += 0.5
                elif c == ' ' or c == '*':
                    bit -= 1
    return bit


if __name__ == '__main__':
    query = args.query.upper().replace('U', 'T')
    p = query_pattern(query, args.errors, args.inserts, args.deletions,
                      args.substitutions, args.wobble)
    if args.mirna:
        p = mirna_pattern(query)
    index = 0
    for f in SeqIO.parse(args.fasta, 'fasta'):
        l = len(str(f.seq))
        for m in p.finditer(str(f.seq).upper()):
            target = str(f.seq[m.start():m.end()]).upper()
            q, t, cig, bit = align(query[::-1], target)
            if args.mirna:
                bit = mirna_score(cig)
                if bit < args.mirna:
                    continue
            energy = 'energy=%.2f;' % fold(query, target) if args.fold else ''
            print('%s\tcoiffR\tmatch\t%d\t%d\t%d\t+\t.\tID=co%d;seq=%s;query=%s;aln=%s;'
                  'score=%.1f;before_end=%d;%s' %
                  (f.name, m.start(), m.end()-1, sum(m.fuzzy_counts),
                   index, t.replace('T', 'U'), q.replace('T', 'U'), cig, bit, (l-m.end()), energy))
            index += 1
        if args.revcomp:
            for m in p.finditer(str(f.seq.reverse_complement())):
                target = str(f.seq.reverse_complement()[m.start():m.end()]).upper()
                q, t, cig, bit = align(query[::-1], target)
                if args.mirna:
                    bit = mirna_score(cig)
                    if bit < args.mirna:
                        continue
                energy = 'energy=%.2f;' % fold(query, target) if args.fold else ''
                print('%s\tcoiffR\tmatch\t%d\t%d\t%d\t-\t.\tID=co%d;seq=%s;query=%s;aln=%s;'
                      'score=%.1f;before_end=%d;%s' %
                      (f.name, l-m.end(), l-m.start()-1, sum(m.fuzzy_counts), index,
                       t.replace('T', 'U'), q.replace('T', 'U'), cig, bit, m.start(), energy))
                index += 1
