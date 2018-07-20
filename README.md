# Chiasm - circRNA identification and characterization suite
Chiasm provides a number of tools for the robust identification of circRNA in RNA-Seq data and helps with the further characterization of detected circRNAs against the genomic background.

The suite evolved as part of my PhD thesis and scientific article (currently in review) investigating circular RNA in honeybees (*Apis mellifera*) and their role in task dependent development of these social insects. The article is a co-operation between the Universities of [Würzburg](https://www.biozentrum.uni-wuerzburg.de/zoo2/personen/wissenschaftler/markus-thamm-persoenliche-seite/) and [Marburg](https://www.uni-marburg.de/de/fb16/ipc/ag-lechner).

We use a combination of `segemehl`+`testrealign` and `BWA`+`CIRI2` for mapping and detection of chiastically mapping read fragments (hence the name "Chiasm"). Identified circRNAs are then characterized with scripts in this repository especially for sequence specific statistics. Each tool can be used on its own, but we provide a shell script pipeline for the aformentioned study as a recipe for an analysis.

## Requirements
In order to recreate our findings based on the same data, or replicate a full analysis for your own data you will need the following programs and environments.

- [Python3](https://www.python.org/): running various scripts
- [Perl5](https://www.perl.org/): running various scripts
- [SRA-tools](https://ncbi.github.io/sra-tools/): downloading `SRA` data and dumping `FASTQ` format files
- [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic): cleaning reads
- [Samtools](http://www.htslib.org/): `SAM` and `BAM` format file manipulations
- [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/): mapping split reads and detection of circRNAs with included `testrealign`
- [BWA](http://bio-bwa.sourceforge.net/): mapping split reads
- [CIRI2](https://sourceforge.net/projects/ciri/): detection of circRNAs in `BWA` mapping data
- [Subread/featureCounts](http://subread.sourceforge.net/): counting reads per feature for normalization
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/): RNAcofold for intron duplex folding energies
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download): intron sequence complementarity analysis and finding homologous regions for miRNA analysis
- [Clustal Omega](http://www.clustal.org/omega/): alignment of homologous regions in the miRNA analysis
- [WebLogo](https://github.com/WebLogo/weblogo): visualizing consensus sequence around splice site

If you want to reproduce only parts of the analysis, pick the software accordingly.

## Running the Pipeline
The shell scripts `pipeline_[0-5]_*.sh` essentially allow a 1:1 replication of our results. Some manual steps might be missing and especially the mapping steps take quite long. So please do not execute them blindly.

- `0_prepare`: downloads the genome data and annotations
- `1_linear`: tries to identify circRNAs from publicly available RNA-Seq data not enriched by RNase R (this was a preliminary screening before we started our study. I suggest you do the same, if your organism of interest is not proven to exhibit circRNAs yet. However, most of these results turned out to be hardly usable without targeted enrichment.)
- `2_circular`: the main mapping and identification part of our pipeline using `segemehl` and `BWA`+`CIRI2` in parallel
- `3_mirna`: analysis of potentially conserved miRNA binding sites within circRNAs
- `4_orthology`: orthological comparison of circRNAs detected by us with already published ones in *Drosophila melanogaster* and *Bombyx mori*
- `5_introns`: characterization of introns flanking circularized exons in comparison to random exons that showed no evidence of circularization
- `6_methylation`: mapping and evaluating WGSBS data of worker bees to determine methylation status of circRNAs

## Scripts in this suite
In the following all standalone scripts that have been written for reuse are explained.

### circquant.py
Quantification of circRNA-Seq mapping results Reads splice information from
segemehls testrealign or circumjunct to quantify circular junctions against
their linear splicing counterparts and outputs a table of all circular
junctions found across multiple samples.

		positional arguments:
		  bedfiles              multiple BED files from containing splice sites
		
		optional arguments:
		  -h, --help            show this help message and exit
		  -c N, --min_circs N   minimum number of found circular junctions in total
		  -e N, --min_experiments N
				                minimum number of samples containing a circular
				                junction
		  -i N, --min_insert N  minimal insert size between junctions (default=100)
		  -I N, --max_insert N  maximal insert size between junctions (default=10000)
		  -a N, --amend N       merge similar junctions and correct by exon boundary
				                within N nt
		  -g FILE, --gff FILE   genome GFF annotation
		  -s, --ignore_strand   ignore junction strand
		  -S, --correct_strand  determine strand info from annotation
		  -n, --normalize       normalize by total library junction reads
		  -f F [F ...], --factors F [F ...]
				                normalize using these factors
		  -H FILE, --host_counts FILE
				                normalize by host gene featureCounts in each library
		  -v, --verbose         print verbose output to STDERR
		
		usage: circquant.py [-h] [-c N] [-e N] [-i N] [-I N] [-a N] [-g FILE] [-s]
				            [-S] [-n | -f F [F ...] | -H FILE] [-v]
				            bedfiles [bedfiles ...]


### circstats.py
Statistics of circRNAs based on genome annotation

		optional arguments:
		  -h, --help            show this help message and exit
		  -g FILE, --gff FILE   genome GFF annotation
		  -c, --correct_strand  correct strand information according to suitable
				                transcript information
		  -v, --verbose         print verbose output to STDERR
		
		usage: circstats.py [-h] [-g FILE] [-c] [-v]


### junctcut.py
Concatenates junctions based on BED format input and genome annotation

		optional arguments:
		  -h, --help            show this help message and exit
		  -v, --verbose         print progress info to STDERR
		  -f FILE, --fasta FILE
				                reference genome sequence
		  -l N, --length N      output length of [exonic] sequence inside junction
				                (default: -1 for full seq)
		  -a N, --amend N       merge similar junctions [and correct by exon boundary]
				                within N nt
		  --print_length_only   output only flanking intron lengths
		
		annotation options:
		  -g FILE, --gff FILE   GFF formated reference annotation
		  -s, --ignore_strand   ignore junction strand info
		  -u N, --upstream N    output 5`-flanking [intronic] sequence of length N[-1
				                for complete intron] (default: 0)
		  -d N, --downstream N  output 3`-flanking [intronic] sequence of length N [-1
				                for complete intron] (default: 0)
		
		usage: junctcut.py [-h] [-v] -f FILE [-l N] [-a N] [-g FILE] [-s] [-u N]
				           [-d N] [--print_length_only]


### circortho.py
Parses ODB relations and annotates a list with orthologs from a second

		positional arguments:
		  orthologs             list of orthologs with ID in first column
		
		optional arguments:
		  -h, --help            show this help message and exit
		  -v, --verbose         print verbose logging info
		  -i INPUT_COLUMN, --input_column INPUT_COLUMN
				                ID field/colum in input
		  -o ORTHO_COLUMN, --ortho_column ORTHO_COLUMN
				                ID field/colum in orthologs
		
		usage: circortho.py [-h] [-v] [-i INPUT_COLUMN] [-o ORTHO_COLUMN] orthologs


### mockgen.py
Generate mock circRNAs based on genome annotation

		optional arguments:
		  -h, --help            show this help message and exit
		  -g FILE, --gff FILE   genome GFF annotation
		  -n N, --number N      number of junctions (default: 100)
		  -i, --internal        exclude first and last exons
		  -s FILE, --similar_length FILE
				                TAB separated file with 5'-/3'-flanking intron length
		  -l N, --length N      minimal length of flanking introns
		  -v, --verbose         print verbose output to STDERR
		
		usage: mockgen.py [-h] -g FILE [-n N] [-i] [-s FILE] [-l N] [-v]


### coiffr.py
Complementary Oligonucleotide Interaction Finder For RNA Searches for regions
in the reference FASTA that are complementary to the query sequence. Allowing
wobble base pairs between G:U, insertions, deletions and substitutions
(errors), if set in the parameters

		optional arguments:
		  -h, --help            show this help message and exit
		  -f FASTA, --fasta FASTA
				                FASTA file of reference genome
		  -q SEQ, --query SEQ   oligonucleotide sequence
		  -r, --revcomp         also search in reverse complement of reference
		  -e N, --errors N      max. allowed errors
		  -i N, --inserts N     max. allowed inserts
		  -d N, --deletions N   max. allowed deletions
		  -s N, --substitutions N
				                max. allowed substitutions
		  -w, --wobble          allow wobble base pairs
		  -m N, --mirna N       miRNA target prediction mode with min. N matches
		  -F, --fold            include free energy score
		  -v, --verbose         print verbose output to STDERR
		
		usage: coiffR [-h] -f FASTA -q SEQ [-r] [-e N] [-i N] [-d N] [-s N] [-w]
				      [-m N] [-F] [-v]

