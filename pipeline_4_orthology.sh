
sed -r 's/^chr(\S+):([0-9]+)-([0-9]+)\s([-\+])\s(\S+)/\1\t\2\t\3\t\5\t\4/' ashwal.tsv | awk 'BEGIN{FS="\t"; OFS="\t"}{$6=$5; split($9, score, ", "); sum=0; for(i=0; i<length(score); i++){sum+=score[i]}; if(sum>=10){$5=sum; $2=$2+1; print}}' | cut -f1-6 > dmel/ashwal.best.bed
circstats -v -g genome/dmel540.gff < dmel/ashwal.best.bed > dmel/ashwal.best.stat.bed

awk 'BEGIN{FS=","; getline; print "#ref\tstart\tend\tname\tscore\tstrand\t5annot\t3annot\tgene\tsymbol"} !/^#/ && $6>=10 {if($4=="-"){t=$7;$7=$8;$8=t}; $2=$2-1; sub("chr", "", $1); print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4"\t"$7"\t"$8"\t"$9"\t"$10}' westholm.csv | sed "s/' UTRs/UTR/g" > dmel/westholm.best.bed
circstats -v -g genome/dmel540.gff < dmel/westholm.best.bed > dmel/westholm.best.stat.bed

# Orthology
wget http://www.orthodb.org/v9/download/odb9_{OG2genes,genes,species}.tab.gz
zgrep mellifera odb9_species.tab.gz     # yields '7460' as the taxid of A.mellifera
zgrep melanogaster odb9_species.tab.gz  # yields '7227' as the taxid of D.melanogaster
zgrep mori odb9_species.tab.gz          # yields '7091' as the taxid of B.mori
zgrep -E '7460:|7227:|7091:' odb9_OG2genes.tab.gz | sort -k2,2 > orthology/odb9_og.tsv     # extract only ortholog groups relevant for amel and dmel
zgrep -E '7460:|7227:|7091:' odb9_genes.tab.gz | sed 's/\t+/\t/g' | sort -k1,1 > orthology/odb9_genes.tsv     # extract only genes relevant for amel dmel, remove empty columns
join -t$'\t' -1 2 -2 1 orthology/odb9_og.tsv orthology/odb9_genes.tsv | sort -k2,2 > orthology/odb9_genes_og.tsv
grep 7460: orthology/odb9_genes_og.tsv | cut -f2,4 > orthology/odb9_genes_og_amel.tsv     # extract only genes relevant for amel
grep 7227: orthology/odb9_genes_og.tsv | cut -f2,4 > orthology/odb9_genes_og_dmel.tsv
grep 7091: orthology/odb9_genes_og.tsv | cut -f2,4 > orthology/odb9_genes_og_bmor.tsv
join -t$'\t' orthology/odb9_genes_og_amel.tsv orthology/odb9_genes_og_dmel.tsv | cut -f2,3 | sort -u | sort -k2,2 > orthology/odb9_amel2dmel_gn.tsv   # join bee with fly genes
join -t$'\t' orthology/odb9_genes_og_amel.tsv orthology/odb9_genes_og_bmor.tsv | cut -f2,3 | sort -u | sort -k2,2 > orthology/odb9_amel2bmor_gn.tsv   # join bee with silk genes
cut -f4,7 dmel/ashwal.best.stat.bed | sed -n '/intragenic/!p' | sort -k2,2 > orthology/ashwal.best.id_gn.tsv 	# ashwal IDs with GeneName
join -t$'\t' -1 2 -2 2 orthology/odb9_amel2dmel_gn.tsv orthology/ashwal.best.id_gn.tsv | cut -f2,3 | sort -u | sort -k1,1 > orthology/ashwal.best.bee2circ.tsv 	# only circRNAs with orthologs
cut -f4,11 dmel/westholm.best.stat.bed | sed -n '/intragenic/!p' | sort -k2,2 > orthology/westholm.best.id_gn.tsv 	# westholm IDs with GeneName
join -t$'\t' -1 2 -2 2 orthology/odb9_amel2dmel_gn.tsv orthology/westholm.best.id_gn.tsv | cut -f2,3 | sort -u | sort -k1,1 > orthology/westholm.best.bee2circ.tsv # only circRNAs with orthologs
sed -rn 's/^.*ID=([^;]+);.*BEEBASE:([^;,]+)[;,].*$/\1\t\2/p' genome/amel.gff | sort -k1,1 -u > genome/amel.id2beebase.tsv
join -t$'\t' -a 1 -2 10 genome/amel.id2beebase.tsv <(sort -k10,10 circ/circ.junctions.bed) | sort -k2,2 > circ/circ.junctions.beebase.tsv 	# enrich circRNAs with BeeBase IDs based on transcript_id
join -t$'\t' -2 2 <(sort -k1,1 -u orthology/odb9_amel2dmel_gn.tsv) circ/circ.junctions.beebase.tsv | cut -f2 | sort -u > circ/circ.junctions.dmel_host_ortho
awk 'BEGIN{OFS="\t"}{if(/^#/){$1="#id\t"$1}else{printf "ame_circ_%8d\t"$0, count++}' circ/circ.junctions.beebase.tsv | sort -k3,3 | join -1 3 - <(sort -k1,1 orthology/odb9_amel2dmel_gn.tsv) | sort -k2,2 -u | wc -l  # 122 circRNA with homolog host genes to fly

circortho -v -i 2 orthology/ashwal.best.bee2circ.tsv < circ/circ.junctions.beebase.tsv > circ/circ.junctions.beebase.ashwal.tsv
circortho -v -i 2 orthology/westholm.best.bee2circ.tsv < circ/circ.junctions.beebase.ashwal.tsv > circ/circ.junctions.beebase.ashwal.westholm.tsv

join -t$'\t' -1 2 -2 2 orthology/odb9_amel2bmor_gn.tsv orthology/gan.id_gn.tsv | cut -f2,3 | sort -u | sort -k1,1 > orthology/gan.bee2circ.tsv 	# only circRNAs with orthologs
sed -rn 's/^.*\sID=([^;]*);.*gene=([^;]*);?.*$/\1\t\2/p' genome/amel.gff | sort -k1,1 -u > genome/amel.id2gene.tsv
join -t$'\t' -a 2 -o auto -e'*' genome/amel.id2gene.tsv <(sort -k1,1 circ/circ.junctions.beebase.ashwal.westholm.tsv) | sort -k2,2 > circ/circ.junctions.beebase.ashwal.westholm.gene.tsv
circortho -v -i 3 orthology/gan.bee2circ.tsv < circ/circ.junctions.beebase.ashwal.westholm.gene.tsv > circ/circ.junctions.beebase.ashwal.westholm.gene.gan.tsv
awk '/\(.*\)/{match($0, /\((\S*)\)/, a); reg="N"; if($5=="Up"){reg="F"}; print a[1]"\t"reg}' liu.dge.tsv > liu.dge.gn2reg.tsv
join -t$'\t' -2 2 -a 2 -o auto -e'*' liu.dge.gn2reg.tsv <(sort -k2,2 circ/circ.junctions.beebase.ashwal.westholm.gene.gan.tsv) > circ/circ.junctions.beebase.ashwal.westholm.gene.gan.liu.tsv

# join IDs into one big list
grep ">" genome/amel.fna | sed -r 's/>(\S*) .*(LG[0-9]*|unplaced|mitochondrion)[ ,].*$/\1\t\2/;s/unplaced/UnGr/;s/mitochondrion/MT/' > genome/amel.chromosome.ids
join -t$'\t' -2 5 -a 2 -o auto -e'*' <(sort -k1,1 genome/amel.chromosome.ids) <(sort -k5,5 circ/circ.junctions.beebase.ashwal.westholm.gene.gan.liu.tsv) > circ.full.tsv
