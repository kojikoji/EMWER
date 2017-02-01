#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -M tiisaishima@gmail.com
#$ -l mem_free=200G
WDIR=/home/tiisaishima/Projects/wfeed/tmp/mimicree
SDIR=/home/tiisaishima/Projects/wfeed/script/mimicree
cd $WDIR
mkdir nsam0${S}_r6_p${P}_g${GN}
cd nsam0${S}_r6_p${P}_g${GN}
mkdir snps
mkdir sim
mkdir sim/sims
for i in {1..12}
do
    python $SDIR/mihap_gen.py 4000 ${P}
    paste $WDIR/hap_ref.mimhap hap${P}.mimhap > hap_init_${P}.mimhap
    python $SDIR/toysoft/pick-random-addtive-snps.py -n 400 -s 0.${S} -e 0.5 -m 0.8 --loci-count 4000 --input hap_init_${P}.mimhap > snps/s0${S}_r6_n$i.txt
    cp  $WDIR/sim_ref.sync sim/s0${S}_r6_n$i.sync
    for j in {12,22,28,45,52,60}
    do
	java -Xmx1g -jar  $SDIR/toysoft/MimicrEESummary.jar --haplotypes-g0 hap_init_${P}.mimhap --recombination-rate ../dmel.rr.txt --output-mode ${j} --replicate-runs 1 --output-format sync --threads 2 --additive snps/s0${S}_r6_n$i.txt --output-file sim/sims/temp_m$j.sync 
    done
    for j in {12,22,28,45,52,60}
    do 
	awk -v GEN=${j} '{print $(NF-1) "\t" $NF "\t" GEN}' sim/sims/temp_m$j.sync >  sim/sims/temp_m$j.dat
	paste sim/s0${S}_r6_n$i.sync sim/sims/temp_m$j.dat > sim/cs0${S}_r6_n$i.sync
	cp sim/cs0${S}_r6_n$i.sync sim/s0${S}_r6_n$i.sync
    done
done
