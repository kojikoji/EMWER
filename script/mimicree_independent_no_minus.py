#!/usr/env python
#-*- coding: utf-8 -*-
import random
import re
import argparse
import os
from cproc.cproc import *
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make sync into dat')
    parser.add_argument('--tflag', '-t',default=False, action = "store_true",help='test flag not qsub')
    parser.add_argument('--num', '-n',default=3,type = int,help='iteration number')
    parser.add_argument('--tsim', '-i',default="tmp/tmp_sim.txt",type = str,help='file name of tmp simulation')
    parser.add_argument('--hap', '-l',default="data/hap_test_one.mimhap",type = str,help='file name of haplotype file')
    parser.add_argument('--recomb', '-c',default="data/dmel.rr.txt",type = str,help='file name of recombination rate file')
    parser.add_argument('--gen', '-g',default="0,30",type = str,help='list of generation')
    parser.add_argument('--rep', '-r',default="1",type = str,help='number of replicate')
    parser.add_argument('--pop', '-p',default="500",type = str,help='number of population')
    parser.add_argument('--bound', '-b',default="0.1,0.9",type = str,help='lower bound of initial allele frequency') # 
    parser.add_argument('--snp', '-v',default="tmp/test/snp_test.snp",type = str,help='list of generation')
    parser.add_argument('--slc', '-s',default="0.1",type = str,help='list of generation')
    parser.add_argument('--afs', '-a',default="",type = str,help='file of AF spectrum')
    parser.add_argument('--output', '-o',type=str,help='haplotype file name')
    args = parser.parse_args()
    cp = cproc(False,args.tflag)
    #make empty file
    if os.path.exists(args.output)==True:
        cp.add("rm " + args.output)
    #make command of haplotype
    mkhapcmd = "python script/make_hap.py"
    mkhapcmd = " ".join([mkhapcmd,"-p",args.pop])
    mkhapcmd = " ".join([mkhapcmd,"-b",args.bound])
    mkhapcmd = " ".join([mkhapcmd,"-o",args.hap])
    mkhapcmd = " ".join([mkhapcmd,"-v",args.snp])
    mkhapcmd = " ".join([mkhapcmd,"-s",args.slc])
    mkhapcmd = " ".join([mkhapcmd,"-a",args.afs])
    #make command of mimicree
    mimcrcmd = "java -Xmx1g -jar  script/mimicree/toysoft/MimicrEESummary.jar --output-format sync --threads 1"
    mimcrcmd = " ".join([mimcrcmd,"--recombination-rate",args.recomb])
    mimcrcmd = " ".join([mimcrcmd,"--haplotypes-g0",args.hap])    
    mimcrcmd = " ".join([mimcrcmd,"--output-mode",args.gen])    
    mimcrcmd = " ".join([mimcrcmd,"--replicate-runs",args.rep])    
    mimcrcmd = " ".join([mimcrcmd,"--additive",args.snp])
    mimcrcmd = " ".join([mimcrcmd,"--output-file",args.tsim])
    #make command of summarization of results
    adrltcmd = "cat"
    adrltcmd = " ".join([adrltcmd,args.tsim,'>>',args.output])    
    for i in range(args.num):
        #Make Haplotype
        cp.add(mkhapcmd)
        #Run simulation 
        cp.add(mimcrcmd)
        #Add result to all results file
        cp.add(adrltcmd)
    cp.exe()
