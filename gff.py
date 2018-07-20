from collections import namedtuple


Feature = namedtuple('Feature', ('chro, type, start, end, strand, sense, name, gene,'
                                 'parents, children, line, attr'))


class GFF:
    filename = ""
    genes = {}
    parents = {}
    features = {}

    def __init__(self, filename, format='auto'):
        self.filename = filename
        with open(filename, 'r') as lines:
            if format == 'auto':
                if (filename.lower().endswith('.gff') or filename.lower().endswith('.gff3')):
                    format = 'gff'
                elif filename.lower().endswith('.gtf'):
                    format = 'gtf'
                else:
                    for l in lines:
                        if 'ID=' in l:
                            format = 'gff'
                            break
                        if 'gene_id ' in l:
                            format = 'gtf'
                            break
                    if format == 'auto':
                        raise Exception('Could not determine GFF/GTF format')
                    lines.seek(0)
            if format.lower() in ['gff', 'gff3']:
                SEP, EQ = ';', '='
                ID, GENE, PARENT = 'id', 'gene', 'parent'
            elif format.lower() == 'gtf':
                SEP, EQ = '; ', ' '
                ID, GENE, PARENT = 'gene_id', 'gene', 'transcript_id'
            else:
                raise Exception('Unknown Genome Annotation Format')

            last_gene = ''
            LN = 0
            for l in lines:
                LN += 1
                if l.startswith('#'):
                    continue
                l = l.rstrip()
                c = l.split('\t')
                if len(c) < 9:
                    continue
                chro, typ, start, end, strand, attr = c[0], c[2], int(c[3]), int(c[4]), c[6], c[8]
                if chro not in self.genes:
                    self.genes[chro] = {'+': [], '-': []}

                attr = {k.lower(): v.strip('"') for k, v in
                        [a.split(EQ, 1) for a in attr.split(SEP) if a]}
                if ID not in attr:
                    continue
                    raise Exception('no ID found for line %d: %s' % (LN, l))
                name = attr[ID]

                gene = attr[GENE] if GENE in attr else ''
                if typ == 'gene':
                    parents = []
                    last_gene = name
                    self.genes[chro][strand].append(name)
                else:
                    parents = attr[PARENT].split(',') if PARENT in attr else [last_gene]
                    for p in parents:
                        self.parents[name] = p
                        if p not in self.features:
                            self.features[p] = Feature(chro, 'parent', 0, 0, '+', True, p,
                                                       gene, 'unknown', [], '', attr)
                        self.features[p].children.append(name)

                children = []
                if name in self.features:
                    children = self.features[name].children
                    if typ == 'CDS':
                        start = min(start, self.features[name].start)
                        end = max(end, self.features[name].end)
                self.features[name] = Feature(chro, typ, start, end, strand, strand == '+',
                                              name, gene, parents, children, l, attr)

    def get_genes(self, chro, strand=None, start=None, end=None):
        g = self.genes[chro]['+']+self.genes[chro]['-'] if strand is None else \
            self.genes[chro][strand]
        g = [self.features[x] for x in g]
        if start is None:
            return g
        if end is None:
            end = start
        return [x for x in g if x.start <= start <= x.end or
                x.start <= end <= x.end or
                start <= x.start <= x.end <= end]

    def get_transcripts_from_gene(self, gene):
        return [self.features[t] for t in gene.children if self.features[t].type == 'mRNA']

    def get_exons_from_transcript(self, transcript, sort_by_direction=True):
        exons = [self.features[x] for x in transcript.children if self.features[x].type == 'exon']
        return sorted(exons, key=lambda x: x.start,
                      reverse=(not transcript.sense and sort_by_direction))

    def get_CDS_from_transcript(self, transcript):
        return [self.features[x] for x in transcript.children if self.features[x].type == 'CDS']

    def get_annotation_from_transcript(self, transcript, pos):
        if pos < transcript.start or transcript.end < pos:
            return 'intragenic'
        if not [e for e in self.get_exons_from_transcript(transcript) if e.start <= pos <= e.end]:
            return 'intronic'
        cds = self.get_CDS_from_transcript(transcript)
        if not cds:
            return 'exonic'
        start = min([c.start for c in cds])
        end = max([c.end for c in cds])
        if pos < start and transcript.sense or not transcript.sense and end < pos:
            return '5UTR'
        if pos < start and not transcript.sense or transcript.sense and end < pos:
            return '3UTR'
        return 'CDS'

    def get_genes_by_feature(self, feature):
        return self.get_genes(feature.chro, feature.start, feature.end, feature.strand)

    def overlap(self, start, end, genes):
        return [self.features[g] for g in genes if start <= self.features[g].start < end or
                start < self.features[g].end <= end or
                self.features[g].start <= start <= end <= self.features[g].end]

    def is_exonic(self, chro, start, end=None, strand=None):
        if end is None:
            end = start
        for g in self.overlap(start, end, self.get_genes(chro, strand)):
            for t in g.children:
                for e in t.children:
                    if e.type == 'exon' and e.start <= start <= e.end:
                        start_found = True
                    if e.type == 'exon' and e.start <= end <= e.end:
                        end_found = True
        return start_found & end_found

    def get_transcripts(self, chro, start, end=None, strand=None):
        if end is None:
            end = start
        # return [f for f in self.features.values() if f.type == 'mRNA' and f.start <= start <= f.end or f.start <= end <= f.end]
        return [self.features[t] for g in self.get_genes(chro, strand, start, end)
                for t in self.features[g.name].children if self.features[t].type in ['mRNA', 'transript']]

# return only exons of one transcript if start and end are part of the transcripts exons
    def get_exons(self, chro, start, end, strand=None):
        best = []
        best_t = None
        for g in self.get_genes(chro, strand, start, end):
            for t in self.features[g.name].children:
                ret, start_exon, end_exon = [], False, False
                for en in self.features[t].children:
                    e = self.features[en]
                    if e.type == 'exon' and (start <= e.start <= end or
                                             start <= e.end <= end):
                        ret.append(e)
                        if e.start <= start <= e.end:
                            start_exon = True
                        if e.start <= end <= e.end:
                            end_exon = True

                if start_exon & end_exon:
                    best = sorted(ret, key=lambda x: x.start)
                    best_t = self.features[t]
                    if start == best[0].start and end == best[-1].end:
                        return (best, self.features[t])
        return (best, best_t)

    def exon_start(self, chro, pos):
        return [e for e in self.get_exons(chro, pos, pos)
                if (e.sense and e.start == pos) or (not e.sense and e.end == pos)]

    def exon_end(self, chro, pos):
        return [e for e in self.get_exons(chro, pos, pos)
                if (not e.sense and e.start == pos) or (e.sense and e.end == pos)]

    def get_left_intron(self, feature):
        last = 0
        exons = [self.features[e] for p in feature.parents for e in self.features[p].children]
        exons = [e for e in exons if e.type == 'exon']
        for e in sorted(exons, key=lambda x: x.start):
            if e.name == feature.name:
                if last == 0:
                    return None
                return feature._replace(type='intron', start=last, end=e.start-1)
            last = e.end + 1

    def get_right_intron(self, feature):
        last = 0
        exons = [self.features[e] for p in feature.parents for e in self.features[p].children]
        exons = [e for e in exons if e.type == 'exon']
        for e in sorted(exons, key=lambda x: x.start, reverse=True):
            if e.name == feature.name:
                if last == 0:
                    return None
                return feature._replace(type='intron', start=e.end+1, end=last)
            last = e.start - 1
