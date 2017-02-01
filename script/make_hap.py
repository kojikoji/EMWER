#!/usr/env python
#-*- coding: utf-8 -*-
import random
import re
import argparse
from util_wfe import *
tag = "2L\t300001\tC\t"
def make_hap(al1,al2,pop,lb=0.1,hb=0.9,afs_file=""):
    su = Spectrum_uniform()
    su.set_low_high(lb,hb)
    if afs_file == "":
        afs_vec = np.ones(50)
    else:
        afs_vec = np.loadtxt(afs_file)
    su.set_spectrum(afs_vec)
    alprob = su.sample(lb,hb)
    alnum = int(alprob*(pop*2))
    als = al1*alnum + al2*(pop*2-alnum)
    als = ''.join(random.sample(als,len(als)))
    als = re.sub(r'(..)',r'\1 ',als)
    als = tag + "C/A\t"+als+"\n"
    return(als)
def make_snp(selection):
    return(tag + selection + '\t' + '0.5')
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make sync into dat')
    parser.add_argument('--population', '-p',type = str,help='population')
    parser.add_argument('--selection', '-s',default = '0.1',type = str,help='selection')
    parser.add_argument('--bound', '-b',type = str,help='boundary of allele frequency: lower,upper')
    parser.add_argument('--output', '-o',type=str,help='haplotype file name')
    parser.add_argument('--snp', '-v',default="",type = str,help='file name of selected snp')
    parser.add_argument('--afs', '-a',default="",type = str,help='file name of allele frequency spectrum')
    args = parser.parse_args()
    with open(args.output,'w') as file:
        lower = float((args.bound).split(',')[0])
        upper = float((args.bound).split(',')[1])
        file.write(make_hap('A','C',int(args.population),lower,upper,args.afs))
    if args.snp != "":
        with open(args.snp,'w') as file:
            file.write(make_snp(args.selection))
