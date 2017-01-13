# -*- coding: utf-8 -*-
from cproc.cproc import cproc
import sys
import os
import csv
if __name__=="__main__":
    sflag = sys.argv[len(sys.argv)-1] == "sim"
    cp = cproc(False)
    #read from condtion file
    f = open("script/test_cond.txt",'rb')
    dread = csv.reader(f,delimiter = '\t')
    dsc = dread.next()[0]
    dnum = dread.next()[0]
    pop = dread.next()[0]
    gen = dread.next()[0]
    slc = dread.next()[0]
    fld = dread.next()[0]
    spnum = dread.next()[0]
    opt_flg = dread.next()[0]
    if os.path.exists("tmp")!=True:
        os.mkdir("tmp")
    fnsim = "tmp/sim_test.dat"
    fnsimc = "tmp/simc_test.dat"
    fnrlt = "tmp/rlt_test.dat"
    fnafs = "script/afs_test.afs"
    cmdsim = "./build/default/WF_sim" \
            + " -n " + spnum \
            + " -r " + dnum \
            + " -p " + pop \
            + " -g " + gen \
            + " -s " + slc \
            + " -f " + fld \
            + " -o " + fnsim \
            + " -c " + fnsimc
    cmdrlt = "./build/default/WF_main" \
            + " -p " + pop \
            + " -d " + dsc \
            + " -s " + "0" \
            + " -i " + fnsim \
            + " -o " + fnrlt \
            + " -q " + fnafs \
            + " -" + opt_flg
    if opt_flg == "l":
        cmdrlt += " -a"
    if sflag :
        cp.add(cmdsim)
    cp.add(cmdrlt)
    cp.exe()
#for debugger
#gdb -i=mi ~/Projects/wfe/bin/WF_test
#r -p 1200 -g 60 -s 0.005 -i ~/Projects/wfe/tmp/in.txt -c ~/Projects/wfe/tmp/inc.txt -o ~/Projects/wfe/tmp/out.txt 
#            + " -c " + fnsimc \
