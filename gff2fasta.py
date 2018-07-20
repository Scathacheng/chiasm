#!/usr/bin/python3

from Bio import SeqIO
import re
import sys
import argparse
from collections import namedtuple


parser = argparse.ArgumentParser('extracts sequences from FASTA based on GFF '
                                 'annotations')
parser.add_argument('-f', '--fasta', help='FASTA genome file', required=True)
parser.add_argument('-g', '--gff', help='GFF annotation file',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-p', '--padding', help='extracted sequences are padded'
                    'by N nucleotides', type=int, default=0, metavar='N')
parser.add_argument('-l', '--max-length', help='limit sequenes to N '
                    'nucleotides', type=int, default=0, metavar='N')
parser.add_argument('-t', '--feature', help='filter by feature in the 3rd '
                    'column', nargs='*', default=[])
parser.add_argument('-i', '--ident', help='identity field in last column',
                    default='ID')
parser.add_argument('-r', '--region', nargs=2, type=int, metavar='START,END',
                    help='extract sequence from START to END or region '
                    '(higher number first for anti-sense strand)')


Seq = namedtuple('Seq', 'chro, start, end, strand, ident, seq')


ID_PATTERN = re.compile('^(\S+)\s+\S+\s+(\S+)\s+'         # chro, _, type
                        '(\d+)\s+(\d+)\s+\S+\s+([\+-])'   # from, to, _, strand
                        '\s+\S+\s+ID=([^;]+);?(.+)$')     # _, ID, rest
GFF_PATTERN = re.compile('^(\S+)\t+[^\t]+\s+(\S+)\s+(\d+)\s+'
                         '(\d+)\s+\S+\s+([\+-]).*$')
GENOME_PATTERN = re.compile('gi\|\d+\|ref\|(NC_\d+\.\d+)\|\s*.*')


def parse_gff_line(line):
    if line.startswith('#'):
        return None
    m = GFF_PATTERN.match(line)
    if m:
        chro = m.group(1)
        feature = m.group(2)
        start = int(m.group(3))
        end = int(m.group(4))
        strand = m.group(5)
        ident = '%s|%d-%d' % (chro, start, end)
        m = re.search(ID_PATTERN, line)
        if m:
            ident = m.group(6)
        if end < start:
            start, end = (end, start)
        return (chro, feature, start, end, strand, ident)


def parse_gff(fasta_name, gff=sys.stdin, padding=0, maxlen=0, typefilter=[]):
    fasta_seqs = SeqIO.to_dict(SeqIO.parse(fasta_name, 'fasta'))

    # sanitizing chromosome names
    tmp = dict()
    changed = 0
    for f in fasta_seqs:
        m = GENOME_PATTERN.match(f)
        if m:
            tmp[m.group(1)] = fasta_seqs[f]
            sys.stderr.write('found GI for refseq %s\n' % f)
            changed = 1
    if changed == 1:
        fasta_seqs = tmp

    if gff == sys.stdin:
        sys.stderr.write("GFF input:")
    for line in gff:
        sys.stderr.write('\r           \r')
        gff_line = parse_gff_line(line)
        if gff_line is None:
            continue
        (chro, feature, start, end, strand, ident) = gff_line
        start -= padding
        end += padding
        if start < 1:
            start = 1
        if end > len(str(fasta_seqs[chro].seq)):
            end = len(str(fasta_seqs[chro].seq))
        if not typefilter or feature in typefilter:
            seq = fasta_seqs[chro].seq[start-1:end]
            if strand is '-':
                seq = seq.reverse_complement()
            if maxlen == 0 or abs(start-end) < maxlen:
                yield Seq(chro, start, end, strand, ident, seq)


def print_fasta(seq):
    print(">%s|%d-%d|%s|%s\n%s\n" % (seq.chro, seq.start, seq.end, seq.strand,
                                     seq.ident, seq.seq))


def main():
    args = parser.parse_args()
    global ID_PATTERN
    ID_PATTERN = re.compile('^(\S+)\s+\S+\s+(\S+)\s+'         # chro, _, type
                            '(\d+)\s+(\d+)\s+\S+\s+([\+-])'   # from, to, _, strand
                            '\s+\S+\s+.*%s=([^;]+);?(.+)$' %  # _, ID, rest
                            args.ident)

    # only a region is specified
    if args.region:
        r = SeqIO.read(args.fasta, 'fasta')
        (start, end, strand) = (args.region[0], args.region[1], '+')
        if start > end:
            seq = str(r.seq[end-1:start].reverse_complement())
            strand = '-'
        else:
            seq = str(r.seq[start-1:end])
        print_fasta(Seq(r.id, start, end, strand, 'region', seq))
        exit()

    for g in parse_gff(args.fasta, args.gff, args.padding, args.max_length,
                       args.feature):
        print_fasta(g)


if __name__ == "__main__":
    main()
