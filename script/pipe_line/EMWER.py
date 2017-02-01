#-*- coding: utf-8 -*-
import sys
sys.path.append("../cproc")
from cproc import cproc 
import os
import argparse
from sim_line_independent import Sim_line_independent
from condition_controler import Condtion_controler
from sync_divider import Sync_divider
from list_mkdir import List_mkdir
from estimator import Estimator
from file_namer import File_namer
class Conducter_sim(Sim_line_independent,Condtion_controler,Sync_divider,List_mkdir,Estimator,File_namer):
    #and root directory of this experiment
    def __init__(self,root_dir='',cp=cproc(False,tflag=True),unit=3):
        self.root_dir = root_dir
        self.cp = cp
        self.unit = unit
        #list of directory keys
        self.list_dir_key = ['target','hap','snp','ref','sim','label','est','afs','dat']
        #list of estimation method, if you want to make comparison with other method, you should see script/estimator.py
        self.est_keys = ['wfe']
        self.est_iter = 100
    #conduct all experiment for condition list
    def conducter(self,cond_file,hap_file,sync_file):
        tag = cond_file.split("/")[-1]
        #make drectory root root/wfe root/sync root/dat etc
        self.list_mkdir()
        #make dictionary and list(key-val_key-val..) of condtion 
        self.make_cond_list(cond_file)
        #loop for condtion
        for cond in self.cond_list :
            #set condtion
            self.set_cond(cond)
            if sync_file == '':
                #simulation is done
                sync_file = self.make_name('sim',cond,'sync')
                self.do_sim(hap_file,sync_file,cond)
            #convert to dat file
            dat_file = self.make_name('dat',cond,'dat')
            self.sync_to_dat(sync_file,dat_file,cond)
            self.est_init(self.est_keys)
            for key in self.est_keys:
                #make cond est
                self.set_est(key,self.cond)
                destf = self.make_name('est',self.cond_est,'est')
                #estimate
                self.do_est(key,dat_file,destf,self.cond_est)
            self.cp.exe()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Simulation for condition and plot roc and AUC')
    parser.add_argument('--tflag', '-t',default=False, action = "store_true",help='test flag not qsub')
    parser.add_argument('--condition', '-c',default='script/cond/emwer_test_slc.cond',type=str,help='condition filename')
    parser.add_argument('--hap_file', '-hp',default='tmp/hap_test_one.mimhap',type=str,help='file of haplotype')
    parser.add_argument('--afs_file', '-a',default='data/afs_unif.afs',type=str,help='file of allele frequency distribution')
    parser.add_argument('--sync_file', '-s',default='',type=str,help='file of allele frequency distribution')
    args = parser.parse_args()
    cp = cproc(False,args.tflag)
    cond_file = args.condition
    root_dir = "tmp/" + (cond_file.split(".")[0]).split("/")[-1] + "/"
    cdt = Conducter_sim(root_dir,cp,unit = 5)
    cdt.conducter(cond_file,args.hap_file,args.sync_file)
    print(cproc(False))
