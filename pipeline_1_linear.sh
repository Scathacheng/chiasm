# FUNCTION: convert SAM to BAM format, sort and index it.
sam2bam() { samtools view -buSh - | samtools sort - > $1; samtools index $1; }

# download and map publicly available SRA data for circRNA screening
cd linear/reads
cat ../../sra_data.tsv | while read sra paired pop tissue enrichment; do
    fastq-dump --gzip --split-files $sra
	t="SE $sra.fastq.gz $sra.qt.fq.gz"
	q="$sra.qt.fq.gz"
    if [[ $paired -eq "paired" ]]; then
        t="PE ${sra}_1.fastq.gz ${sra}_2.fastq.gz ${sra}_1.qt.fq.gz ${sra}_1.unpaired.fq.gz ${sra}_2.qt.fq.gz ${sra}_2.unpaired.fq.gz"
        q="${sra}_1.qt.fq.gz -p ${sra}_2.qt.fq.gz"
    fi
	trimmomatic $t SLIDINGWINDOW:3:15 MINLEN:20
    segemehl -i ../../genome/amel.mehl.idx -d ../../genome/amel.fna -q $q -s -S -t 12 -E 0.1 -w 1.0 | sam2bam $sra.bam
	testrealign -d ../../genome/amel.fna -n -q <(samtools view -h $sra.bam) --transfile ../$sra.trans.bed --splitfile ../$sra.split.bed
    rm ${sra}*.gz
done
cd -
featureCounts -T 8 -a genome/amel.gtf -o linear/reads/featureCounts.tsv linear/reads/*.bam
circquant -v -a 5 -i 50 -I 50000 -S -g genome/amel.gff -H linear/reads/featureCounts.tsv -s linear/*.split.bed > linear/linear.quant.bed
circstats -c -v -g genome/amel.gff < linear/linear.quant.bed | tee linear/linear.quant.stat.bed | wc -l # 6266 junctions (6214 with exon boundary)
awk '/^#/ || $5 >= 10{print}' linear/linear.quant.stat.bed | tee linear/linear.quant10.stat.bed | wc -l	# 381 with min 10 reads
cp linear/linear.quant10.stat.bed linear/linear.junctions.bed
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u 10 -d 10 < linear/linear.junctions.bed > linear/linear.junctions.splice.fna # 372 exonic
grep -i ........AGGT linear/linear.junctions.splice.fna | wc -l 	# 359 with usual splice signal
seq2logo < linear/linear.junctions.splice.fna > linear/linear.junctions.splice.ps

# Intron analysis
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u 500 -d 500 < linear/linear.junctions.bed > linear/linear.intron500.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 5 -u 500 -d 500 < linear/linear.junctions.bed > linear/linear.intron505.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u -1 < linear/linear.junctions.bed > linear/linear.5intron.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -d -1 < linear/linear.junctions.bed > linear/linear.3intron.fna
awk '/^>/ {getline; print length($0)}' < linear/linear.5intron.fna > linear/linear.5intron.lengths
awk '/^>/ {getline; print length($0)}' < linear/linear.3intron.fna > linear/linear.3intron.lengths


