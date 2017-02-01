# EMWER
##Introduction
EM algorithm of Wright-Fisher model for analyzing Evolve and Resequence data
##Instllation
Please, check configuration using the following command,
and install the required libraries.

```shell-session
$ ./waf configure
```

Please, install using the following command

```shell-session
$ ./waf build
```

The executable file will be generated as `build/default/WF_main`

##Usage
You can execute our program through python script `script/EMWER.py`.
By default, simulation of ER data by mimicree and
estimation of selection coefficients are performed.

```shell-session
$ python script/EMWER.py
```

The estimation result is written in

```
tmp/emwer_test_slc/est/opt-slc_bd-0.05,0.95_num-2_dep-30_rep-3_pop-250_disc-30_afs-unif_chrom-2L_gen-0,30_slc-0.05_est-wfe.est
```

The sample is like below:

```
250     0.0439589       2.84145 33
250     -0.0518147      2.81102 40
```

Col1,2,3,4 respectively represent effective poplation size, selection coefficient, likelihood ratio, estimation time[ms].

The experimental condition is drawn in `script/cond/emwer_test_slc.cond` like below:

```
gen:0,30
pop:250
slc:0.05
dep:30
rep:3
disc:30
opt:slc
num:2
bd:0.05,0.95
afs:unif
```

Each row represents each experimental condition.;`gen` is the generation of sampling points in experimental evolution, `pop` is effective population size, `slc` is selection coefficients in simulation, `dep` is mean of sequence depth, `rep` is the number of population replicates, `dis` is the number of discrete states of allele frequecy used in estimation, `opt` is the name of parameter which is optimized, `num` is the number of SNPs contained in ER data, `bd` is the upper and lower bound of initial allele frequency, `afs` is the type of allele frequency spectorum used in simulation. If you want to  experiment on multiple conditions at once, you can write tab separated multiple values for each row like `gen:0,30\t0,60`.

You can test the estimation of effective population size by:

```shell-session
$ python script/EMWER.py -c script/cond/emwer_test_pop.cond
```

`-c` specify the file name of experimental conndition.

Furthermore,you can test the estimation
```shell-session
$ python script/EMWER.py -c script/cond/emwer_test_real.cond -s data/real_data_test.sync
```
`-s` specify the synchronized file (see detail in https://sourceforge.net/p/popoolation2/wiki/Tutorial/) of ER data[Orozco-terWengel,2012].

The experimental condition is little bit changed like:

```
pop:200
repl:1,4,5,4,5,1,5,4,5,1
genl:0,0,0,15,15,23,27,37,37,37
slc:0
opt:slc
disc:30
rbd:0.05,0.95
afs:default
```

`genl` represents the list of generation when corresponding colmn in the synchronized file was sampled, and `repl` represents the list of replicates index of the colmn.
