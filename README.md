# Chiasm - circRNA identification and characterization suite
Chiasm provides a number of tools for the robust identification of circRNA in RNA-Seq data and helps with the further characterization of detected circRNAs against the genomic background.

The suite evolved as part of my PhD thesis and scientific article (currently in review) investigating circular RNA in honeybees (*Apis mellifera*) and their role in task dependent development of these social insects. The article is a co-operation between the Universities of [Würzburg](https://www.biozentrum.uni-wuerzburg.de/zoo2/personen/wissenschaftler/markus-thamm-persoenliche-seite/) and [Marburg](https://www.uni-marburg.de/de/fb16/ipc/ag-lechner).

We use a combination of `segemehl`+`testrealign` and `BWA`+`CIRI2` for mapping and detection of chiastically mapping read fragments (hence the name "Chiasm"). Identified circRNAs are then characterized with scripts in this repository especially for sequence specific statistics. Each tool can be used on its own, but we provide a shell script pipeline for the aformentioned study as a recipe for an analysis.

## Requirements
In order to recreate our findings based on the same data, or replicate a full analysis for your own data you will need the following programs and environments.

- Python3
- Perl5
- sratools
- trimmomatic
- samtools
- segemehl
- BWA
- CIRI2
- subread/featureCounts
- ViennaRNA
- BLAST

I you want to reproduce only parts of the analysis, pick the software accordingly.

## Running the Pipeline
The shell scripts `pipeline_[0-5]_*.sh` essentially allow a 1:1 replication of our results. Some manual steps might be missing and especially the mapping steps take quite long. So please do not execute them blindly.

- `0_prepare`: downloads the genome data and annotations
- `1_linear`: tries to identify circRNAs from publicly available RNA-Seq data not enriched by RNase R (this was a preliminary screening before we started our study. I suggest you do the same, if your organism of interest is not proven to exhibit circRNAs yet. However, most of these results turned out to be hardly usable without targeted enrichment.)
- `2_circular`: the main mapping and identification part of our pipeline using `segemehl` and `BWA`+`CIRI2` in parallel
- `3_mirna`: analysis of potentially conserved miRNA binding sites within circRNAs
- `4_orthology`: orthological comparison of circRNAs detected by us with already published ones in *Drosophila melanogaster* and *Bombyx mori*
- `5_introns`: characterization of introns flanking circularized exons in comparison to random exons that showed no evidence of circularization.
