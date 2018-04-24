# EMWER
## Introduction
This is the repository for an EM algorithm of Wright-Fisher model for analyzing Evolve and Resequence (ER) data.
From ER data, which contain allele frequency paths during experimental evolution for numerous SNPs, this program can estimate population genetics parameters i.e.) effective population size, selection coefficients and dominance parameters for each SNP or all SNPs.

## Installation
Please, check configuration using the following command,
and install the required libraries.

```shell-session:
$ ./waf configure
```

Please, install using the following command

```shell-session
$ ./waf build
```

The executable file will be generated as `build/default/WF_main`

## Usage
Although you can execute the above executable file directly, we recommend executing our program through python script `EMWER.py` as below:

```shell-session
$ python EMWER.py -c script/cond/emwer_test_real.cond -s data/real_data_test.sync
```
The `-s` specifies the `synchronized file` (see detail in https://sourceforge.net/p/popoolation2/wiki/Tutorial/) of ER data. The `-c` specifies the file name of the `condition file`. The `-o` specifies the file name of the `estimation result`.


## Estimation result
The `estimation result` for two SNPs is like below:

```
250	0.0364121	0.5	0.000328661	0	3.84428	-1
250	0.0130407	0.5	0.000610262	0	0.278235	-1
```

Col1,2,3,4,5,6,7 respectively represent effective population size, selection coefficient, dominance parameter, half width of 95% confidence interval for  estimated selection coefficient, that for dominance parameter, likelihood ratio, original sign of selection coefficients for variant allele. We noticed that, as default, selection coefficients and dominance parameter is those of the allele which is positively selected. Hence, if you want to derive these parameter for variant allele, you have to specify the condition `allele` as `variant` in `condition file` discussed below.


## Condition for estimation
Experimental and estimation conditions can be specified in `condition file` as below:

```
genl:0,0,0,15,15,23,27,37,37,37
repl:1,4,5,4,5,1,5,4,5,1
pop:200
slc:0
dom:0.5
opt:s
disc:200
rbd:0.05,0.95
allele:positive
```


Each row can specify each experimental condition.;

`genl` specifies the list of generation when each column in the synchronized file was sampled.

`repl` specifies the list of replicate index of each column in the synchronized.

`pop`, `slc` and `dom` specifies initial values of estimation for effective population size, selection coefficient and dominance parameter, respectively.

`disc` specifies the maximum number of discrete states of allele frequency used in estimation.

`opt` specifies the first letter of parameter name which you want to optimize (`p`, `s` and `d` means effective population sizes, selection coefficients and dominance parameter respectively; you can conduct simultaneous estimation for selection and dominance by specifying `opt` as `sd`).

`rbd` specifies the upper and lower limit of allele frequency of population at 0th generation.

`allele` specifies the allele which selection coefficients and dominance parameters are defined for. When you specify `allele` as `positive`, EMWER calculate the parameters for positively selected alleles, which results in positive values for estimated selection coefficients. When you specify `allele` as `variant`, EMWER calculate the parameters for variant alleles.



