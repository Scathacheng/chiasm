# 1. RNAcofold 5' vs 3'

# circ
perl ./Cofold_energy.pl circ/circ.5intron.fna circ/circ.3intron.fna 50
perl ./Cofold_energy.pl circ/circ.5intron.fna circ/circ.3intron.fna 100
perl ./Cofold_energy.pl circ/circ.5intron.fna circ/circ.3intron.fna 200

RNAfold -C --noPS <cofoldS.50.tmp | tee cofoldS.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldS.50.mfe
RNAfold -C --noPS <cofoldSE.50.tmp | tee cofoldSE.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldSE.50.mfe
RNAfold -C --noPS <cofoldES.50.tmp | tee cofoldES.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldES.50.mfe
RNAfold -C --noPS <cofoldE.50.tmp | tee cofoldE.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldE.50.mfe

RNAfold -C --noPS <cofoldS.100.tmp | tee cofoldS.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldS.100.mfe
RNAfold -C --noPS <cofoldSE.100.tmp | tee cofoldSE.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldSE.100.mfe
RNAfold -C --noPS <cofoldES.100.tmp | tee cofoldES.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldES.100.mfe
RNAfold -C --noPS <cofoldE.100.tmp | tee cofoldE.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldE.100.mfe

RNAfold -C --noPS <cofoldS.200.tmp | tee cofoldS.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldS.200.mfe
RNAfold -C --noPS <cofoldSE.200.tmp | tee cofoldSE.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldSE.200.mfe
RNAfold -C --noPS <cofoldES.200.tmp | tee cofoldES.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldES.200.mfe
RNAfold -C --noPS <cofoldE.200.tmp | tee cofoldE.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldE.200.mfe
for f in cofold*.{out,mfe}; do mv $f circ.$f; done

# mock
perl ./Cofold_energy.pl mock/mock.5intron.fna mock/mock.3intron.fna 50
perl ./Cofold_energy.pl mock/mock.5intron.fna mock/mock.3intron.fna 100
perl ./Cofold_energy.pl mock/mock.5intron.fna mock/mock.3intron.fna 200

RNAfold -C --noPS <cofoldS.50.tmp | tee cofoldS.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldS.50.mfe
RNAfold -C --noPS <cofoldSE.50.tmp | tee cofoldSE.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldSE.50.mfe
RNAfold -C --noPS <cofoldES.50.tmp | tee cofoldES.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldES.50.mfe
RNAfold -C --noPS <cofoldE.50.tmp | tee cofoldE.50.out | perl parseMFE.pl | sort -n -k1,1 > cofoldE.50.mfe

RNAfold -C --noPS <cofoldS.100.tmp | tee cofoldS.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldS.100.mfe
RNAfold -C --noPS <cofoldSE.100.tmp | tee cofoldSE.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldSE.100.mfe
RNAfold -C --noPS <cofoldES.100.tmp | tee cofoldES.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldES.100.mfe
RNAfold -C --noPS <cofoldE.100.tmp | tee cofoldE.100.out | perl parseMFE.pl | sort -n -k1,1 > cofoldE.100.mfe

RNAfold -C --noPS <cofoldS.200.tmp | tee cofoldS.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldS.200.mfe
RNAfold -C --noPS <cofoldSE.200.tmp | tee cofoldSE.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldSE.200.mfe
RNAfold -C --noPS <cofoldES.200.tmp | tee cofoldES.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldES.200.mfe
RNAfold -C --noPS <cofoldE.200.tmp | tee cofoldE.200.out | perl parseMFE.pl | sort -n -k1,1 > cofoldE.200.mfe
for f in cofold*.{out,mfe}; do mv $f mock.$f; done



### OLD!!!
#   vvvvvv

threaded.pl Cofold_energy.job
for F in cofold?.out; do perl ./parseMFE.pl <$F | sort -n >$F.circ.mfe ; done
for F in cofold?.out; do perl ./parseMFE.pl <$F | sort -n | awk '{print NR " " $0}' >$F.circ.mfe ; done

# mock
perl ./Cofold_energy.pl mock/mock.5intron.fna mock/mock.3intron.fna
threaded.pl Cofold_energy.job
for F in cofold?.out; do perl ./parseMFE.pl <$F | sort -n >$F.mock.mfe ; done
for F in cofold?.out; do perl ./parseMFE.pl <$F | sort -n | awk '{print NR " " $0}' >$F.mock.mfe ; done

R

cA = read.table("cofoldA.out.circ.mfe")
cB = read.table("cofoldB.out.circ.mfe")
cC = read.table("cofoldC.out.circ.mfe")
cD = read.table("cofoldD.out.circ.mfe")

mA = read.table("cofoldA.out.mock.mfe")
mB = read.table("cofoldB.out.mock.mfe")
mC = read.table("cofoldC.out.mock.mfe")
mD = read.table("cofoldD.out.mock.mfe")

boxplot(c(cA,cB,cC,cD,mA,mB,mC,mD))

plot(c(cA,cB,cC,cD,mA,mB,mC,mD))


# BLAST bit score of RCMs in the first and last 100nt
perl blast_score.pl circ/circ.5intron.fna circ/circ.3intron.fna
perl blast_score.pl mock/mock.5intron.fna mock/mock.3intron.fna

circ.blast.starts = read.table(paste0(basedir, 'circ.blast_S.tmp'), header = F, stringsAsFactors = F)
circ.blast.ends = read.table(paste0(basedir, 'circ.blast_E.tmp'), header = F, stringsAsFactors = F)
circ.blast.se = read.table(paste0(basedir, 'circ.blast_SE.tmp'), header = F, stringsAsFactors = F)
circ.blast.es = read.table(paste0(basedir, 'circ.blast_ES.tmp'), header = F, stringsAsFactors = F)
mock.blast.starts = read.table(paste0(basedir, 'mock.blast_S.tmp'), header = F, stringsAsFactors = F)
mock.blast.ends = read.table(paste0(basedir, 'mock.blast_E.tmp'), header = F, stringsAsFactors = F)
mock.blast.se = read.table(paste0(basedir, 'mock.blast_SE.tmp'), header = F, stringsAsFactors = F)
mock.blast.es = read.table(paste0(basedir, 'mock.blast_ES.tmp'), header = F, stringsAsFactors = F)
boxplot(c(circ.blast.starts, mock.blast.starts, circ.blast.ends, mock.blast.ends, circ.blast.se, mock.blast.se, circ.blast.es, mock.blast.es),
		at=c(1,2, 5,6, 9,10, 13,14))
