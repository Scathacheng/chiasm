# miRNA analysis
wget ftp://mirbase.org/pub/mirbase/21/mature.fa.gz
zcat mature.fa.gz | sed -n '/mellifera/{p;n;p}' > mirna/amel.mirna.fna
cp mature.fa.gz mirna/mature.fa.gz

junctcut -v -f genome/amel.fna -g genome/amel.gff -l -1 < circ/circ.junctions.bed > circ/circ.junctions.fna
sed -nr 's/^>\S+\s(id.*)$/\1/;s/ /\n/gp' circ/circ.junctions.fna > circ/circ.exons.ids
sort -u circ/circ.exons.ids | while read id; do grep "ID=$id;" genome/amel.gff; done | gff2fasta -f genome/amel.fna > circ/circ.exons.fna
sed 's/^>//;s/ .*$//;N;s/\n/\t/' mirna/amel.mirna.fna | while read name mirna; do coiffr -f circ/circ.exons.fna -q $mirna -F -m 15 > mirna/circ.$name.coiffr.gff; done   # use coiffr in miRNA mode (score mostly starting positions)
for f in mirna/circ.*.coiffr.gff; do if [ $(cat $f | wc -l) == '0' ]; then rm $f; fi; done  # delete empty files
for f in mirna/circ.*.coiffr.gff; do mirna=$(echo $f | sed 's/^.*\/circ.//;s/\..*$//'); sed 's/match/'$mirna'/' $f; done > mirna/circ.mirna.targets.gff   # concatenate results into one file
awk 'BEGIN{FS="\t"; OFS="\t"} $7=="+" {match($1, "^(.*)\\|([0-9]+)-([0-9]+)\\|([-//+])\\|(.*)$", a); $1=a[1]; $2=a[5]; if(a[4]=="+") {$4+=a[2]; $5+=a[2]} else {if($7=="+") {$7="-"} else {$7="+"}; start=$4; $4=a[3]-$5; $5=a[3]-start}; print}' mirna/circ.mirna.targets.gff > mirna/circ.mirna.gff

junctcut -v -f genome/amel.fna -g genome/amel.gff -l -1 < mock/mock.junctions.bed > mock/mock.junctions.fna
sed -nr 's/^>\S+\s(id.*)$/\1/;s/ /\n/gp' mock/mock.junctions.fna > mock/mock.exons.ids
sort -u mock/mock.exons.ids | while read id; do grep "ID=$id;" genome/amel.gff; done | gff2fasta -f genome/amel.fna > mock/mock.exons.fna
sed 's/^>//;s/ .*$//;N;s/\n/\t/' mirna/amel.mirna.fna | while read name mirna; do coiffr -f mock/mock.exons.fna -q $mirna -F -m 15 -r > mirna/mock.$name.coiffr.gff; done   # use coiffr in miRNA mode (score mostly starting positions)
for f in mirna/mock.*.coiffr.gff; do if [ $(cat $f | wc -l) == '0' ]; then rm $f; fi; done  # delete empty files
for f in mirna/mock.*.coiffr.gff; do mirna=$(echo $f | sed 's/^.*\/mock.//;s/\..*$//'); sed 's/match/'$mirna'/' $f; done > mirna/mock.mirna.targets.gff   # concatenate results into one file
awk 'BEGIN{FS="\t"; OFS="\t"} $7=="+" {match($1, "^(.*)\\|([0-9]+)-([0-9]+)\\|([-//+])\\|(.*)$", a); $1=a[1]; $2=a[5]; if(a[4]=="+") {$4+=a[2]; $5+=a[2]} else {if($7=="+") {$7="-"} else {$7="+"}; start=$4; $4=a[3]-$5; $5=a[3]-start}; print}' mirna/mock.mirna.targets.gff > mirna/mock.mirna.gff

# Conserved miRNA targets

#1. +/-1 100 nt ->  fastacmd
awk 'BEGIN{OFS="\t"} {$4+=-100; $5+=100; print $0;}' mirna/circ.targets.gff > mirna/circ.targets.extended.gff
./fastacmd-pro.mod.pl -type=gff genome/amel.fna mirna/circ.targets.extended.gff > mirna/circ.targets.extended.fasta

#2. Find conserved regions in other genomes with BLAST
for F in amel ador aflo acer edil lven mqua bimp bter dmel bmor; do
    makeblastdb -dbtype nucl -in genome/$F.fna;
    blastn -num_threads 16 -db genome/$F.fna -query mirna/circ.targets.extended.fasta -evalue 1e-05 -outfmt 6 | sort -k12,12nr | ./bestonly.pl > mirna/circ.targets.$F.bla;
    fastacmd-pro.pl -type=blast genome/$F.fna mirna/circ.targets.$F.bla > mirna/circ.targets.$F.fna;
done

mkdir -p mirna/target_groups/
cd mirna/target_groups/

# Group sequences by same target into individual files
rm *.fasta; for F in ../circ.targets.????.fna; do perl ../fna2group.pl $F <$F ; done

# Alignment of files with more than one sequence (single sequences are AMEL only)
for F in *.fasta; do [[ $(wc -l $F | cut -d" " -f1) -lt 4 ]] || clustalo -i $F >$F.aln; done
cd -

# Gather all alignments and identify which hits are conserved among apis, eusocial insects or insects in general
for F in mirna/target_groups/*.aln; do perl getTargetSite.pl $F circ/circ.quant.stat.id.5fold.bed; done > mirna/circ.conserved.targets.tsv


### SOME STATS
grep -v ">" ../circ/circ.exons.fna | wc -c				#> 96413 nucleotides in circRNA bounding exons
cat circ.exons.seed6.genomic.gff | sort -k1,1 -k4,5 -u | wc -l  # number of unique miRNA binding sites
awk '$7=="A"{print}' circ.conserved.tsv | sort -k1,3 -u | wc -l # number of unique conserved Apis sites
awk '$8=="E"{print}' circ.conserved.tsv | sort -k1,3 -u | wc -l # number of unique conserved Eusocial sites

# control
perl ../mock_exons.pl < ../mock/mock.exons.fna;
rm mock.exons.43.fna # delete last, because it is only half full
for i in $(seq 1 42); do cat conserved/amel.mirna.seed6.tsv | while read name mirna; do coiffr -f mock.exons.$i.fna -q $mirna -F | sed 's/match/'$name'/'; done > mock.exons.$i.targets; done
for i in $(seq 1 42); do awk 'BEGIN{FS="\t"; OFS="\t"} $7=="+" {match($1, "^(.*)\\|([0-9]+)-([0-9]+)\\|([-//+])\\|(.*)$", a); $1=a[1]; $2=a[5]; if(a[4]=="+") {$4+=a[2]-100; $5+=a[2]+100} else {if($7=="+") {$7="-"} else {$7="+"}; $4=a[3]-$5-100; $5=a[3]-$4+100}; print}' mock.exons.$i.targets > mock.targets.$i.genomic.extended.gff; ./fastacmd-pro.mod.pl -type=gff ../genome/amel.fna mock.targets.$i.genomic.extended.gff > mock.targets.$i.fna; done
for i in $(seq 1 42); do
	for F in conserved/????.fna; do g=$(echo $F | sed 's/.*\///;s/.fna//'); echo $g; blastn -num_threads 16 -db $F -query mock.targets.$i.fna -evalue 1e-05 -outfmt 6 | sort -k12,12nr | perl bestonly.pl > mock.$i.$g.bla ; done
	for G in conserved/????.fna; do g=$(echo $G | sed 's/.*\///;s/.fna//'); echo $g; fastacmd-pro.pl -type=blast $G mock.$i.$g.bla > mock.$i.$g.bla.fna; done
	rm group.*.fasta; for F in mock.$i.????.bla.fna; do perl fna2group.pl $F; done
	rm group.*.aln; for F in group.*.fasta; do [[ $(wc -l $F | cut -d" " -f1) -lt 4 ]] || clustalo -i $F >$F.aln; done
	for F in group.*.fasta.aln; do perl getTargetSite.pl $F <(echo ""); done > mock.$i.conserved.tsv
done
rm group.*.fasta
rm group.*.aln


