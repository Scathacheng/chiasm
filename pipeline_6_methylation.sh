# Methylation analysis of Apis mellifera genomes

## Getting the Data
#| Run       | Sample       | Assay         | Source  | Mbases  | Mbytes  |
#|-----------|--------------|---------------|---------|---------|---------|
#| SRR445809 | WorkerPool5  | Bisulfite-Seq | GENOMIC | 8,183   | 5,782   |
#| SRR445808 | WorkerPool4  | Bisulfite-Seq | GENOMIC | 7,152   | 5,014   |
#| SRR445807 | WorkerPool3  | Bisulfite-Seq | GENOMIC | 10,559  | 6,944   |
#| SRR445806 | WorkerPool2  | Bisulfite-Seq | GENOMIC | 6,553   | 4,277   |
#| SRR445805 | WorkerPool1  | Bisulfite-Seq | GENOMIC | 5,602   | 3,689   |
#| SRR445799 | NursePool1   | Bisulfite-Seq | GENOMIC | 24,139  | 14,654  |
#| SRR445778 | NursePool2   | Bisulfite-Seq | GENOMIC | 14,082  | 9,926   |
#| SRR445777 | NursePool3   | Bisulfite-Seq | GENOMIC | 15,342  | 10,741  |
#| SRR445776 | NursePool4   | Bisulfite-Seq | GENOMIC | 24,317  | 17,152  |
#| SRR445775 | NursePool5   | Bisulfite-Seq | GENOMIC | 26,406  | 18,462  |
#| SRR445774 | NursePool6   | Bisulfite-Seq | GENOMIC | 6,079   | 4,253   |
#| SRR445773 | ForagerPool1 | Bisulfite-Seq | GENOMIC | 15,639  | 10,185  |
#| SRR445771 | ForagerPool2 | Bisulfite-Seq | GENOMIC | 17,026  | 11,083  |
#| SRR445770 | ForagerPool3 | Bisulfite-Seq | GENOMIC | 15,720  | 10,981  |
#| SRR445769 | ForagerPool4 | Bisulfite-Seq | GENOMIC | 7,377   | 5,276   |
#| SRR445768 | ForagerPool5 | Bisulfite-Seq | GENOMIC | 12,286  | 8,732   |
#| SRR445767 | ForagerPool6 | Bisulfite-Seq | GENOMIC | 8,953   | 6,270   |
cd methylation
cat data.tsv | while read sra rest; do 
	fastq-dump --gzip --split-files $sra;
done
cd -

# Simplified genome only based on continuous chromosomes (some tools will crash otherwise)
mkdir genome/bismark
awk '/>NW_/{exit}{print}' genome/amel.fna > genome/bismark/amel.chr.fna


## Bismark
bismark_genome_preparation --bowtie2 genome/bismark
cat methylation/data.tsv | while read sra rest; do
	bismark -p 14 --bowtie2 -n 1 genome/bismark -1 methylation/${sra}_1.fastq.gz -2 methylation/${sra}_2.fastq.gz
	bismark_methylation_extractor -s --gzip --cytosine_report --bedGraph --genome_folder genome/bismark --comprehensive ${sra}_1_bismark_bt2_pe.bam
done


## Evaluation
# Extract methylation info only for circRNAs in the HIGH set with separate 5' and 3'-introns
tail -n+2 circ/circ.junctions.bed | cut -f1-6,29,30 | while read chro s e name circ strand sintron eintron; do echo ">$name"; for f in methylation/*.cov.gz; do echo ">$name $chro $s $strand $e $circ" >> $f.exon; zcat $f | awk '$1=="'$chro'" && '$s'<$2 && $2<'$e' {print}' >> $f.exon; done; done
tail -n+2 circ/circ.junctions.bed | cut -f1-6,29,30 | while read chro start end name circ strand sintron eintron; do e=$start; s=$(($start-$sintron)); if [ "$strand" == "-" ]; then s=$end; e=$(($end+$sintron)); fi; echo ">$name $start $end $sintron $eintron $strand $s = $e"; for f in methylation/*.cov.gz; do echo ">$name $chro $s $strand $e $circ" >> $f.5intron; zcat $f | awk '{if($1=="'$chro'" && '$s'<$2 && $2<'$e') {print; passed=1} else if(passed) {exit}}' >> $f.5intron; done; done
tail -n+2 circ/circ.junctions.bed | cut -f1-6,29,30 | while read chro start end name circ strand sintron eintron; do s=$end; e=$(($end+$eintron)); if [ "$strand" == "-" ]; then s=$(($start-$eintron)); e=$start; fi; echo ">$name $start $end $sintron $eintron $strand $s = $e"; for f in methylation/*.cov.gz; do echo ">$name $chro $s $strand $e $circ" >> $f.3intron; zcat $f | awk '{if($1=="'$chro'" && '$s'<$2 && $2<'$e') {print; passed=1} else if(passed) {exit}}' >> $f.3intron; done; done


# Generate methylation stats for different lengths around exon
./methyl_stats.py -c 5 -m 0.1 -v methylation/SR*.exon > methylation/circ.exon
./methyl_stats.py -c 5 -m 0.1 -v methylation/SR*.5intron > methylation/circ.5intron
./methyl_stats.py -c 5 -m 0.1 -v methylation/SR*.3intron > methylation/circ.3intron
./methyl_stats.py -f 50 -c 5 -m 0.1 -v methylation/SR*.3intron > methylation/circ.50.3intron
./methyl_stats.py -l 50 -c 5 -m 0.1 -v methylation/SR*.5intron > methylation/circ.50.5intron
./methyl_stats.py -f 100 -c 5 -m 0.1 -v methylation/SR*.3intron > methylation/circ.100.3intron
./methyl_stats.py -l 100 -c 5 -m 0.1 -v methylation/SR*.5intron > methylation/circ.100.5intron
./methyl_stats.py -f 200 -c 5 -m 0.1 -v methylation/SR*.3intron > methylation/circ.200.3intron
./methyl_stats.py -l 200 -c 5 -m 0.1 -v methylation/SR*.5intron > methylation/circ.200.5intron


# Methylation stats for the generated random sample of "linear" exons:
tail -n+2 mock/mock.junctions.bed | cut -f1-6,12,13 | while read chro start end name circ strand sintron eintron; do name="$chro#$start$strand$end"; s=$end; e=$end; echo ">$name $start $end $sintron $eintron $strand $s = $e"; for f in methylation/*.cov.gz; do echo ">$name $chro $s $strand $e $circ" >> $f.mock.exon; zcat $f | awk '{if($1=="'$chro'" && '$s'<$2 && $2<'$e') {print; passed=1} else if(passed) {exit}}' >> $f.mock.exon; done; done
tail -n+2 mock/mock.junctions.bed | cut -f1-6,12,13 | while read chro start end name circ strand sintron eintron; do name="$chro#$start$strand$end"; e=$start; s=$(($start-$sintron)); if [ "$strand" == "-" ]; then s=$end; e=$(($end+$sintron)); fi; echo ">$name $start $end $sintron $eintron $strand $s = $e"; for f in methylation/*.cov.gz; do echo ">$name $chro $s $strand $e $circ" >> $f.mock.5intron; zcat $f | awk '{if($1=="'$chro'" && '$s'<$2 && $2<'$e') {print; passed=1} else if(passed) {exit}}' >> $f.mock.5intron; done; done
tail -n+2 mock/mock.junctions.bed | cut -f1-6,12,13 | while read chro start end name circ strand sintron eintron; do name="$chro#$start$strand$end"; s=$end; e=$(($end+$eintron)); if [ "$strand" == "-" ]; then s=$(($start-$eintron)); e=$start; fi; echo ">$name $start $end $sintron $eintron $strand $s = $e"; for f in methylation/*.cov.gz; do echo ">$name $chro $s $strand $e $circ" >> $f.mock.3intron; zcat $f | awk '{if($1=="'$chro'" && '$s'<$2 && $2<'$e') {print; passed=1} else if(passed) {exit}}' >> $f.mock.3intron; done; done

# Generate methylation stats for different lengths around exon
./methyl_stats.py -c 5 -m 0.1 -v methylation/SR*.mock.exon > methylation/mock.exon
./methyl_stats.py -c 5 -m 0.1 -v methylation/SR*.mock.5intron > methylation/mock.5intron
./methyl_stats.py -c 5 -m 0.1 -v methylation/SR*.mock.3intron > methylation/mock.3intron
./methyl_stats.py -f 50 -c 5 -m 0.1 -v methylation/SR*.mock.3intron > methylation/mock.50.3intron
./methyl_stats.py -l 50 -c 5 -m 0.1 -v methylation/SR*.mock.5intron > methylation/mock.50.5intron
./methyl_stats.py -f 100 -c 5 -m 0.1 -v methylation/SR*.mock.3intron > methylation/mock.100.3intron
./methyl_stats.py -l 100 -c 5 -m 0.1 -v methylation/SR*.mock.5intron > methylation/mock.100.5intron
./methyl_stats.py -f 200 -c 5 -m 0.1 -v methylation/SR*.mock.3intron > methylation/mock.200.3intron
./methyl_stats.py -l 200 -c 5 -m 0.1 -v methylation/SR*.mock.5intron > methylation/mock.200.5intron


# Convert to BED format
for f in methylation/*.exon; do sed -r 's/^>(\S+) (\S+) (\S+) \S+ (\S+) \S+ \S+ \S+ \S+ (\S+)$/>\1 \2 \3 \4 \5/' $f; done


# Create BED file for methylation sites
for f in methylation/*.cov.gz; do grep -vh ">" $f.{exon,5intron,3intron} > $f.circ.bed; done


#Create SEG file from BED to visualize in IGV
for f in methylation/*.circ.bed; do sra=$(echo $f | sed 's/_.*$//'); awk '{print "'$sra'\t"$1"\t"$2"\t"$3"\t"$5"\t"$4}' $f; done > methylation/all.seg

