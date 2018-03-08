# FUNCTION: convert SAM to BAM format, sort and index it.
sam2bam() { samtools view -buSh - | samtools sort - > $1; samtools index $1; }

# mapping circRNA enriched data
cd circ/reads
for sra in SRR4343845 SRR4343846; do # single end
    fastq-dump --gzip --split-files $sra
    segemehl -i ../../genome/amel.mehl.idx -d ../../genome/amel.fna -q $sra.fq.gz -S -s -t 12 -E 0.1 -w 1.0 | sam2bam ../mapping/$sra.bam
    testrealign -d ../../genome/amel.fna -n -q <(samtools view -h ../mapping/$sra.bam) --transfile ../mapping/$sra.trans.bed --splitfile ../mapping/$sra.split.bed
done

for sra in SRR4343847 SRR4343848; do # paired end
    fastq-dump --gzip --split-files $sra
    segemehl -i ../../genome/amel.mehl.idx -d ../../genome/amel.fna -q ${sra}_1.fastq.gz -p ${sra}_2.fastq.gz -S -s -t 12 -E 0.1 -w 1.0 | sam2bam ../mapping/$sra.bam
    testrealign -d ../../genome/amel.fna -n -q <(samtools view -h ../mapping/$sra.bam) --transfile ../mapping/$sra.trans.bed --splitfile ../mapping/$sra.split.bed
done
cd -
featureCounts -T 8 -a genome/amel.gtf -o circ/featureCounts.tsv circ/mapping/*.bam

circquant -v -a 5 -c 10 -i 50 -I 50000 -S -g genome/amel.gff -H circ/featureCounts.tsv -s circ/mapping/*.split.bed > circ/circ.quant.bed
circstats -c -v -g genome/amel.gff < circ/circ.quant.bed | tee circ/circ.quant.stat.bed | wc -l     # 3511 (3487 with exon boundary)
awk 'BEGIN{OFS="\t"}!/^#/{split($4, a, ":"); $7=a[3]":"a[4]; $4=sprintf("ame_circ_%07d", count++)}{print}' circ/circ.quant.stat.bed > circ/circ.quant.stat.id.bed
awk '/^#/ || $5 >= 10{print}' circ/circ.quant.stat.id.bed | tee circ/circ.quant10.stat.bed | wc -l 	# 486 with 10+ reads
awk '/^#/ || $16>5*$13 && $14>2 && ($17>2 || $20>2) {print}' circ/circ.quant10.stat.bed | tee circ/circ.quant.5fold.stat.bed | wc -l 	# 253 enriched 5fold
awk '/^#/ || $16>10*$13 && $14>2 && ($17>2 || $20>2) {print}' circ/circ.quant10.stat.bed | tee circ/circ.quant.10fold.stat.bed | wc -l 	# 196 enriched 10fold
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u 10 -d 10 < circ/circ.quant.5fold.stat.bed > circ/circ.quant.5fold.splice.fna  # 1 without annotation
grep -i ........AGGT circ/circ.quant.5fold.splice.fna | wc -l 	# 247 canonical splice signal
seq2logo < circ/circ.quant.5fold.splice.fna > circ/circ.quant.5fold.splice.ps
cp circ/circ.quant.5fold.stat.bed circ/circ.junctions.bed
join <(cut -f1-3 circ/circ.junctions.bed | sed 's/\t//g' | sort) <(cut -f1-3 linear/circles/linear.quant.bed | sed 's/\t//g' | sort) | wc -l    # 210 also in linear data


junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u 500 -d 500 < circ/circ.junctions.bed > circ/circ.intron500.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 5 -u 500 -d 500 < circ/circ.junctions.bed > circ/circ.intron505.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u -500 < circ/circ.junctions.bed > circ/circ.5intron500.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -d -500 < circ/circ.junctions.bed > circ/circ.3intron500.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -u -1 < circ/circ.junctions.bed > circ/circ.5intron.fna
junctcut -v -g genome/amel.gff -f genome/amel.fna -l 0 -d -1 < circ/circ.junctions.bed > circ/circ.3intron.fna
awk '/^>/ {getline; print length($0)}' < circ/circ.5intron.fna > circ/circ.5intron.lengths
awk '/^>/ {getline; print length($0)}' < circ/circ.3intron.fna > circ/circ.3intron.lengths
