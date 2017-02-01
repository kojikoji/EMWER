#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -M tiisaishima@gmail.com
#$ -l mem_free=200G
WDIR=/home/tiisaishima/Projects/wfeed/tmp/mimicree
HDIR=/home/tiisaishima/Projects/wfeed
SDIR=/home/tiisaishima/Projects/wfeed/script/mimicree
cd $HDIR
make -j4
cd $WDIR
cd nsam0${S}_r6_p${P}_g${GN}
mkdir pdata
for i in {1..12}
do
    $SDIR/sync_process2.sh sim/s0${S}_r6_n$i.sync > pdata/s0${S}_r6_n$i.txt
done
wait
for i in {1..12}
do
    $HDIR/bin/WF_mimicre ${P} ${GN} pdata/s0${S}_r6_n$i.txt pdata/es0${S}_r6_n$i.est &
done
wait
mkdir pred
for i in {1..12}
do
    awk '{print $(NF - 1)-$(NF - 2)}'  pdata/es0${S}_r6_n$i.est >  pred/lhs0${S}_r6_n$i.pred
    if test $i -eq 1
    then
	  cp pred/lhs0${S}_r6_n$i.pred lhs0${S}_r6_p${P}_g${GN}.pred
    else
	  paste clhs0${S}_r6_p${P}_g${GN}.pred pred/lhs0${S}_r6_n$i.pred > lhs0${S}_r6_p${P}_g${GN}.pred
    fi
    head -n 10 lhs0${S}_r6_p${P}_g${GN}.pred
    cp lhs0${S}_r6_p${P}_g${GN}.pred clhs0${S}_r6_p${P}_g${GN}.pred
done
mv lhs0${S}_r6_p${P}_g${GN}.pred simres/		 
