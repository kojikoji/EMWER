#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -M tiisaishima@gmail.com
#$ -l mem_free=200G
WDIR=/home/tiisaishima/Projects/wfeed/tmp/mimicree
cd $WDIR
mkdir cmh
for i in {1..12}; do perl ../toysoft/popoolation2/cmh-test.pl --min-count 1 --min-coverage 1 --max-coverage 100000  --population 1-2,3-4,5-6 --input sim/s0${S}_r6_n$i.sync --output cmh/s0${S}_r6_n$i.cmh; done 
mkdir simres
python ../toysoft/rocr-generate-labellist.py sim/s0${S}_r6_n1.sync snps/s0${S}_r6_n{1..12}.txt > simres/lbs0${S}_r6_p${P}_g${GN}.pred
python ../toysoft/rocr-generate-predictionlist.py cmh/s0${S}_r6_n{1..12}.cmh > simres/cmhs0${S}_r6_p${P}_g${GN}.pred
R --vanilla --args simres/s0${S}.labels simres/s0${S}.predictions simres/s0${S}.labels simres/s0${S}.predictions roc < ../toysoft/toyroc.R
