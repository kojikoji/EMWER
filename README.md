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
You can execute our program through python script `script/pipe_line/EMWER.py`.
By default, simulation of ER data by mimicree and
estimation of selection coefficients are performed.

```shell-session
$ python script/pipe_line/EMWER.py
```

The estimation result is written in

```
tmp/emwer_test_slc/est/opt-slc_bd-0.05,0.95_num-2_dep-30_rep-3_pop-250_disc-100_afs-unif_chrom-2L_gen-0,30_slc-0.05_est-wfe.est
```

The sample is like below:

```
250	0.0364121	0.5	0.000328661	0	3.84428	-1
250	0.0130407	0.5	0.000610262	0	0.278235	-1
```

Col1,2,3,4,5,6,7 respectively represent effective poplation size, selection coefficient, dominance parameter, half width of 95% confidence interval for  estimated selection coeffcient, that for dominace parameter, likelihood ratio, original sign of selection coefficients for variant allele.

The experimental condition is drawn in `script/cond/emwer_test_slc.cond` like below:

```
gen:0,30
pop:250
slc:0.05
dep:30
rep:3
disc:100
opt:s
num:2
bd:0.05,0.95
afs:unif
```

Each row represents each experimental condition.;`gen` is the generation of sampling points in experimental evolution, `pop` is effective population size, `slc` is selection coefficients in simulation, `dep` is mean of sequence depth, `rep` is the number of population replicates, `dis` is the number of discrete states of allele frequecy used in estimation, `opt` is the initial of parameter name which is optimized ("s", "p" and "d" means selection coefficients, effective population sizes and dominance parameter respectively; when you can specify simultaneous estimation for selection and dominance by "sd"), `num` is the number of SNPs contained in ER data, `bd` is the upper and lower bound of initial allele frequency, `afs` is the type of allele frequency spectorum used in simulation. If you want to  experiment on multiple conditions at once, you can write tab separated multiple values for each row like `gen:0,30\t0,60`.

You can test the estimation of effective population size by:

```shell-session
$ python script/pipe_line/EMWER.py -c script/cond/emwer_test_pop.cond
```

`-c` specify the file name of experimental conndition.

Furthermore,you can test the estimation
```shell-session
$ python EMWER.py -c script/cond/emwer_test_real.cond -s data/real_data_test.sync
```
`-s` specify the synchronized file (see detail in https://sourceforge.net/p/popoolation2/wiki/Tutorial/) of ER data[Orozco-terWengel,2012].

The experimental condition is little bit changed like:

```
pop:200
repl:1,4,5,4,5,1,5,4,5,1
genl:0,0,0,15,15,23,27,37,37,37
slc:0
opt:s
disc:100
rbd:0.05,0.95
afs:default
```

`genl` represents the list of generation when corresponding colmn in the synchronized file was sampled, and `repl` represents the list of replicates index of the colmn.
