# EMWER
##Introduction
EM algorithm of Wright-Fisher model for analyzing Evolve and Resequence data
##Instllation
Please, configure using the following command:

```shell-session
$ ./waf configure
```

Please, install using the following command

```shell-session
$ ./waf build
```

The executable file will be generated as `build/default/WF_main`

##Usage
You can execute our program, like bellow 
```shell-session
$ ./build/default/WF_main -p 200 -d 30 -s 0 -i tmp/test/sim_test.dat -o tmp/test/rlt_test.est -q tmp/test/afs_test.afs -f
```
Out put file is `tmp/test/rlt_test.est`
