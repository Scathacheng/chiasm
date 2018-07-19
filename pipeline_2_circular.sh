# FUNCTION: convert SAM to BAM format, sort and index it.
sam2bam() { samtools view -buSh - | samtools sort - > $1; samtools index $1; }


# SEGEMEHL PIPELINE

# generate SEGEMEHL index of the genome (3.7GB)
segemehl -x genome/amel.mehl.idx -d genome/amel.fna

# mapping circRNA enriched data
cd circ/reads
for sra in SRR4343845 SRR4343846; do # single end
    fastq-dump --gzip $sra
    segemehl -i ../../genome/amel.mehl.idx -d ../../genome/amel.fna -q $sra.fastq.gz -S -s -t 12 -E 0.1 -w 1.0 | sam2bam ../mapping/$sra.bam
    testrealign -d ../../genome/amel.fna -n -q <(samtools view -h ../mapping/$sra.bam) --transfile ../mapping/$sra.trans.bed --splitfile ../mapping/$sra.split.bed
done

for sra in SRR4343847 SRR4343848; do # paired end
    fastq-dump --gzip --split-files $sra
    segemehl -i ../../genome/amel.mehl.idx -d ../../genome/amel.fna -q ${sra}_1.fastq.gz -p ${sra}_2.fastq.gz -S -s -t 12 -E 0.1 -w 1.0 | sam2bam ../mapping/$sra.bam
    testrealign -d ../../genome/amel.fna -n -q <(samtools view -h ../mapping/$sra.bam) --transfile ../mapping/$sra.trans.bed --splitfile ../mapping/$sra.split.bed
done
cd -

# counting reads per feature with the `featureCounts` software from the `subread` package
featureCounts -T 8 -a genome/amel.gtf -o circ/featureCounts.tsv circ/mapping/*.bam

# aggregating BSJs from different libraries, joining almost identical junctions within 5nt and justifying them to annotated exon boundaries within 5nt
circquant -v -a 5 -c 1 -i 50 -I 50000 -S -g genome/amel.gff -H circ/featureCounts.tsv -s circ/mapping/*.split.bed > circ/circ.mehl.quant.bed
	# parsed 29005 individual junctions
	# merged 636 junctions ...
	# read 98646 exons (132215 redundant) ...
	# removed 5206 junctions (insert: 5206, few reads: 0, few samples: 0)
	# exonic junction start: 5037 (337 amended) / exonic junction ends: 4994 (317 amended)
	# corrected strand information for 1845 junctions

# generate additional annotation for statistical analysis
circstats -c -v -g genome/amel.gff < circ/circ.quant.bed | tee circ/circ.quant.stat.bed | wc -l     # 3511 (3487 with exon boundary)


# CIRI2 PIPELINE

# mapping with BWA
bwa index -a bwtsw genome/amel.fna
bwa mem -T 19 -t 14 genome/amel.fna circ/reads/SRR4343847_1.fq.gz rcirc/reads/SRR4343847_2.fq.gz > circ/mapping/SRR4343847.bwa.sam
bwa mem -T 19 -t 14 genome/amel.fna circ/reads/SRR4343848_1.fq.gz rcirc/reads/SRR4343848_2.fq.gz > circ/mapping/SRR4343848.bwa.sam
bwa mem -T 19 -t 14 genome/amel.fna circ/reads/SRR4343845.fq.gz > circ/mapping/SRR4343845.bwa.sam
bwa mem -T 19 -t 14 genome/amel.fna circ/reads/SRR4343846.fq.gz > circ/mapping/SRR4343846.bwa.sam

# analysis with CIRI2
for sra in SRR4343847 SRR4343848 SRR4343845 SRR4343846; do
	perl CIRI2.pl -I circ/mapping/$sra.bwa.sam -O circ/mapping/$sra.bwa.ciri2.tsv -F genome/amel.fna -A genome/amel.gtf -T 14
	# Convert from CIRI2 to `testrealign`-like BED:
	awk 'NR==1{$2="#"$2}{print $2"\t"$3"\t"$4"\tsplits:"$5":"($5+$7)":"($5+$7)":C:P\t0\t"$11}' circ/mapping/$sra.bwa.ciri2.tsv > circ/mapping/$sra.bwa.ciri2.bed
done

# aggregate results with `circquant` to make analysis comparable (we do not lose many junctions):
circquant -v -a 5 -c 1 -i 50 -I 50000 -S -g genome/amel.gff -H circ/featureCounts.tsv -s circ/mapping/*.bwa.ciri2.bed > circ/circ.ciri2.quant.bed
	# parsed 1962 individual junctions
	# merged 13 junctions ...
	# read 98646 exons (132215 redundant) ...
	# removed 29 junctions (insert: 29, few reads: 0, few samples: 0)
	# exonic junction start: 1738 (11 amended) / exonic junction ends: 1721 (15 amended)
	# corrected strand information for 823 junctions


# Further Analysis
# The two tables were joined and filtered to only contain BSJs detected by both approaches with at least 3 reads, minimum 10 reads across all libraries in either method and 5-fold enrichment in RNase R treated samples vs untreated
# The result was saved in `BED`-like format in `circ/circ.junctions.bed`

# splice site
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u 10 -d 10 < circ/circ.junctions.bed > circ/circ.junctions.splice.fna  # 1 without annotation
grep -i ........AGGT circ/circ.junctions.splice.fna | wc -l 	# 247 canonical splice signal
seq2logo < circ/circ.quant.5fold.splice.fna > circ/circ.quant.5fold.splice.ps

join <(cut -f1-3 circ/circ.junctions.bed | sed 's/\t//g' | sort) <(cut -f1-3 linear/circles/linear.quant.bed | sed 's/\t//g' | sort) | wc -l    # 210 also in linear data

# extracting flanking introns at different lengths
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u 500 -d 500 < circ/circ.junctions.bed > circ/circ.intron500.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 5 -u 500 -d 500 < circ/circ.junctions.bed > circ/circ.intron505.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u -500 < circ/circ.junctions.bed > circ/circ.5intron500.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -d -500 < circ/circ.junctions.bed > circ/circ.3intron500.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u -1 < circ/circ.junctions.bed > circ/circ.5intron.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -d -1 < circ/circ.junctions.bed > circ/circ.3intron.fna
awk '/^>/ {getline; print length($0)}' < circ/circ.5intron.fna > circ/circ.5intron.lengths
awk '/^>/ {getline; print length($0)}' < circ/circ.3intron.fna > circ/circ.3intron.lengths
