# required software
# sratools
# samtools
# trimmomatic
# segemehl
# subread/featureCounts

# download genomes for further annotation and genomic analysis
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Apis_mellifera/latest_assembly_versions/GCF_000002195.4_Amel_4.5/GCF_000002195.4_Amel_4.5_genomic.fna.gz
zcat GCF_000002195.4_Amel_4.5_genomic.fna.gz > genome/amel.fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Apis_mellifera/latest_assembly_versions/GCF_000002195.4_Amel_4.5/GCF_000002195.4_Amel_4.5_genomic.gff.gz
zcat GCF_000002195.4_Amel_4.5_genomic.gff.gz > genome/amel.gff
wget ftp://ftp.flybase.net/genomes/dmel/dmel_r5.40_FB2011_08/fasta/dmel-all-chromosome-r5.40.fasta.gz
zcat dmel-all-chromosome-r5.40.fasta.gz > genome/dmel540.fna
wget ftp://ftp.flybase.net/genomes/dmel/dmel_r5.40_FB2011_08/gff/dmel-all-r5.40.gff.gz
zcat dmel-all-r5.40.gff.gz > genome/dmel540.gff

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/442/555/GCF_001442555.1_ACSNU-2.0/GCF_001442555.1_ACSNU-2.0_genomic.gff.gz
zcat GCF_001442555.1_ACSNU-2.0_genomic.gff.gz > genome/acer.gff
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/442/555/GCF_001442555.1_ACSNU-2.0/GCF_001442555.1_ACSNU-2.0_genomic.fna.gz
zcat GCF_001442555.1_ACSNU-2.0_genomic.fna.gz > genome/acer.fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/785/GCF_000184785.2_Aflo_1.0/GCF_000184785.2_Aflo_1.0_genomic.gff.gz
zcat GCF_000184785.2_Aflo_1.0_genomic.gff.gz > genome/aflo.gff
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/785/GCF_000184785.2_Aflo_1.0/GCF_000184785.2_Aflo_1.0_genomic.fna.gz
zcat GCF_000184785.2_Aflo_1.0_genomic.fna.gz > genome/aflo.fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/469/605/GCF_000469605.1_Apis_dorsata_1.3/GCF_000469605.1_Apis_dorsata_1.3_genomic.gff.gz
zcat GCF_000469605.1_Apis_dorsata_1.3_genomic.gff.gz > genome/ador.gff
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/469/605/GCF_000469605.1_Apis_dorsata_1.3/GCF_000469605.1_Apis_dorsata_1.3_genomic.fna.gz
zcat GCF_000469605.1_Apis_dorsata_1.3_genomic.fna.gz > genome/ador.fna

# delete compressed genome files
rm *.gz

# convert GFF annotation to GTF conform format for transcriptomic read counting
sed -rn 's/([A-Za-z0-9_-]+)=([^;]+);?\s?/\1 "\2"; /gp' genome/amel.gff > genome/amel.pseudo.gtf
awk '{ if(match($0, "ID \"([^\"]+)\";?", id)) { if($3 == "gene") { last_gene = id[1]; last_transcript = ""}; if($3 == "transcript" || $3 == "mRNA") {last_transcript = id[1]}; for(i=1;i<=NF;i++){if(i==9) printf("gene_id \"%s\"; transcript_id \"%s\"; ", last_gene, last_transcript); printf($i); if(i<9) printf("\t"); if(i>=9) printf(" ")} print ""}}' genome/amel.pseudo.gtf > genome/amel.gtf

# generate SEGEMEHL index of the genome (3.7GB)
segemehl -x genome/amel.mehl.idx -d genome/amel.fna
