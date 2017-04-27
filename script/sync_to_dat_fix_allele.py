#-*- coding: utf-8 -*-
import argparse
import numpy as np
from scipy.stats import binom    
def get_val_list(line,st,ed): 
    list_hap = line.strip().split("\t")
    list_hap = list_hap[st:ed]
    #print(st)
    #print(ed)
    #print(list_hap)
    return(list_hap)
def get_ref(line):
    refl = "\t".join(get_val_list(line,0,3))
    return(refl)
def make_af(syncl,genl):
    sumfrq = np.zeros(len(syncl[0].split(":")))
    i = 0
    for sync in syncl:
        if genl[i] == '0':
            frq = np.array(map(int,sync .split(":")))
            sumfrq += frq
            i += 1
    #primary allele and second allele
    pal = -1
    sal = -1
    pfrq = -1
    sfrq = -1
    sorted_sumfrq = sorted(sumfrq,reverse=True)
    sumfrql = []
    for frq in sumfrq:
        sumfrql.append(frq)
    pal =sumfrql.index(sorted_sumfrq[0])
    sal =sumfrql.index(sorted_sumfrq[1])
    afl = []
    #print(pal)
    #print(sal)
    for sync in syncl:
        frql = sync.split(":")
        afl.append([frql[0],frql[2]])
    return(afl)
# rep-gen が正しくした配列での各要素の元の配列での順番
def ind_list_maker(genl,repl):
    if len(genl) != len(repl) :
        print("ERROR!!:: genl and repl is different length")
        exit(-1)
    gri_lst = []
    for ind in range(0,len(genl)):
        gri = {}
        gri['ind'] = ind
        gri['gen'] = genl[ind]
        gri['rep'] = repl[ind]
        gri_lst.append(gri)
    #先にgen
    sorted_gri_lst = sorted(gri_lst,key = lambda x : x['gen'])
    #次にrepo
    sorted_gri_lst = sorted(sorted_gri_lst,key = lambda x : x['rep'])
    sorted_ind_lst = []
    for gri in sorted_gri_lst:
        sorted_ind_lst.append(gri['ind'])
    return(sorted_ind_lst)
def renew_afs(disc,pfrq,sfrq,old_afs):
    h = float(1)/(disc-1)
    dep = pfrq + sfrq
    ad_afs =  np.zeros(disc)
    sum_lh = 0
    for i in range(disc):
        likelihood = binom.pmf(pfrq,dep,i*h)
        ad_afs[i] = likelihood
        sum_lh += likelihood
    ad_afs /= sum_lh
    return(old_afs + ad_afs)
def sync_to_dat(sync_file,dat_file,ref_file,afs_file,genl,repl,disc,bd):
    lower = float(bd.split(",")[0])
    upper = float(bd.split(",")[1])
    aflag = False
    if afs_file != "":
        aflag = True
    rflag = False
    if ref_file != "":
        rflag = True
        reff = open(ref_file,"w")
    #open file
    syncf = open(sync_file)
    datf = open(dat_file,"w")
    afs = np.zeros(disc)
    afs_count = 0
    #one liner
    for line in syncf:
        #record initial allele frequency and how many 
        iaf = 0
        iafc = 0
        #bdflag become false when initi af is out of bound
        bdflag = True
        oline =""
        prep = '-1'
        afl = make_af(get_val_list(line,3,3+len(genl)),genl)
        #print(afl)
        #print(len(afl))
        ind_list =  ind_list_maker(genl,repl)
        for ind in ind_list:
            #if population is changed tab is inserted
            if repl[ind] != prep:
                oline += "\t"
            else:
                oline += ","
            prep = repl[ind]
            af = afl[ind]
            afgen = af + [genl[ind]]
            oline += ",".join(afgen)
            # sum afs if generation = 0
            if genl[ind] == "0" and aflag:
                pfrq = int(af[0])
                sfrq = int(af[1])
                if max(pfrq,sfrq) < 1:
                    pfrq = 1
                af = float(pfrq)/(pfrq+sfrq)
                if af > upper:
                    bdflag = False
                if af < lower:
                    bdflag = False    
                afs = renew_afs(disc,pfrq,sfrq,afs)
                afs_count += 1
                # for recording initial AF
                iaf += af
                iafc += 1
        oline = oline.strip() + '\n'
        if bdflag :
            datf.write(oline)
            if rflag :
                if iafc != 0:                    
                    iaf = iaf/iafc
                reff.write(get_ref(line) + "\t" + str(iaf) + '\n')
    datf.close()
    syncf.close()
    if aflag:
        afs = afs/afs_count
        np.savetxt(afs_file,afs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make sync into dat')
    parser.add_argument('--genlist', '-g',type = str,help='list of generation sync')
    parser.add_argument('--replist', '-r',type=str,help='list of replication of sync')
    parser.add_argument('--input', '-i',type=str,help='sync file name')
    parser.add_argument('--output', '-o',type=str,help='dat file name')
    parser.add_argument('--afs', '-a',type=str,default="",help='allele frequency distribution')
    parser.add_argument('--ref', '-f',type=str,default="",help='reference file of snps')
    parser.add_argument('--bound', '-b',type=str,default="0,1",help='bouder of allele frequency')
    parser.add_argument('--disc', '-d',type=int,default=50,help='discretization of allele frequency distribution')
    args = parser.parse_args()
    sync_to_dat(args.input,args.output,args.ref,args.afs,args.genlist.split(','),args.replist.split(','),args.disc,args.bound)
