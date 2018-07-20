# BLAST bit score of RCMs in the first and last 100nt
./blast_score.pl circ/circ.5intron.fna circ/circ.3intron.fna
./blast_score.pl mock/mock.5intron.fna mock/mock.3intron.fna


# 1. RNAcofold 5' vs 3'
# circ
./cofold_energy.pl circ/circ.5intron.fna circ/circ.3intron.fna 50
./cofold_energy.pl circ/circ.5intron.fna circ/circ.3intron.fna 100
./cofold_energy.pl circ/circ.5intron.fna circ/circ.3intron.fna 200

RNAfold -C --noPS <cofoldS.50.tmp | tee cofoldS.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldS.50.mfe
RNAfold -C --noPS <cofoldSE.50.tmp | tee cofoldSE.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldSE.50.mfe
RNAfold -C --noPS <cofoldES.50.tmp | tee cofoldES.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldES.50.mfe
RNAfold -C --noPS <cofoldE.50.tmp | tee cofoldE.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldE.50.mfe

RNAfold -C --noPS <cofoldS.100.tmp | tee cofoldS.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldS.100.mfe
RNAfold -C --noPS <cofoldSE.100.tmp | tee cofoldSE.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldSE.100.mfe
RNAfold -C --noPS <cofoldES.100.tmp | tee cofoldES.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldES.100.mfe
RNAfold -C --noPS <cofoldE.100.tmp | tee cofoldE.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldE.100.mfe

RNAfold -C --noPS <cofoldS.200.tmp | tee cofoldS.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldS.200.mfe
RNAfold -C --noPS <cofoldSE.200.tmp | tee cofoldSE.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldSE.200.mfe
RNAfold -C --noPS <cofoldES.200.tmp | tee cofoldES.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldES.200.mfe
RNAfold -C --noPS <cofoldE.200.tmp | tee cofoldE.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldE.200.mfe
for f in cofold*.{out,mfe}; do mv $f circ.$f; done

# random control introns
./cofold_energy.pl mock/mock.5intron.fna mock/mock.3intron.fna 50
./cofold_energy.pl mock/mock.5intron.fna mock/mock.3intron.fna 100
./cofold_energy.pl mock/mock.5intron.fna mock/mock.3intron.fna 200

RNAfold -C --noPS <cofoldS.50.tmp | tee cofoldS.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldS.50.mfe
RNAfold -C --noPS <cofoldSE.50.tmp | tee cofoldSE.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldSE.50.mfe
RNAfold -C --noPS <cofoldES.50.tmp | tee cofoldES.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldES.50.mfe
RNAfold -C --noPS <cofoldE.50.tmp | tee cofoldE.50.out | ./parseMFE.pl | sort -n -k1,1 > cofoldE.50.mfe

RNAfold -C --noPS <cofoldS.100.tmp | tee cofoldS.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldS.100.mfe
RNAfold -C --noPS <cofoldSE.100.tmp | tee cofoldSE.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldSE.100.mfe
RNAfold -C --noPS <cofoldES.100.tmp | tee cofoldES.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldES.100.mfe
RNAfold -C --noPS <cofoldE.100.tmp | tee cofoldE.100.out | ./parseMFE.pl | sort -n -k1,1 > cofoldE.100.mfe

RNAfold -C --noPS <cofoldS.200.tmp | tee cofoldS.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldS.200.mfe
RNAfold -C --noPS <cofoldSE.200.tmp | tee cofoldSE.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldSE.200.mfe
RNAfold -C --noPS <cofoldES.200.tmp | tee cofoldES.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldES.200.mfe
RNAfold -C --noPS <cofoldE.200.tmp | tee cofoldE.200.out | ./parseMFE.pl | sort -n -k1,1 > cofoldE.200.mfe
for f in cofold*.{out,mfe}; do mv $f mock.$f; done

