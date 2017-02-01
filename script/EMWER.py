#-*- coding: utf-8 -*-
import sys
from cproc import cproc 
import os
import argparse
from pipe_line.sim_line_independent import Sim_line_independent
from pipe_line.condition_controler import Condtion_controler
from pipe_line.sync_divider import Sync_divider
from pipe_line.list_mkdir import List_mkdir
from pipe_line.estimator import Estimator
from pipe_line.file_namer import File_namer
class Conducter_sim(Sim_line_independent,Condtion_controler,Sync_divider,List_mkdir,Estimator,File_namer):
    #and root directory of this experiment
    def __init__(self,root_dir='',cp=cproc(False,tflag=True),unit=3):
        self.root_dir = root_dir
        self.cp = cp
        self.unit = unit
        #list of directory keys
        self.list_dir_key = ['target','hap','snp','ref','sim','label','est''afs']
        #list of estimation method, if you want to make comparison with other method, you should see script/estimator.py
        self.est_keys = ['wfe']
        self.est_iter = 100
    #conduct all experiment for condition list
    def conducter(self,cond_file,hap_file,recomb_file):
        tag = cond_file.split("/")[-1]
        #make drectory root root/wfe root/sync root/dat etc
        self.list_mkdir()
        #make dictionary and list(key-val_key-val..) of condtion 
        self.make_cond_list(cond_file)
        self.recomb_file = recomb_file
        #loop for condtion
        for cond in self.cond_list :
            #set condtion
            self.set_cond(cond)
            condhap_file = self.make_name('hap',cond,'mimhap')
            #make interval and chrom select for haplotype
            self.preprocess(hap_file,condhap_file,cond)
            #simulation is done
            sync_file = self.make_name('sim',cond,'sync')
            self.do_sim(condhap_file,sync_file,cond)
            #convert to dat file
            dat_file = self.make_name('dat',cond,'dat')
            self.sync_to_dat(sync_file,dat_file,cond)
            self.est_init(self.est_keys)
            for key in self.est_keys:
                #make cond est
                self.set_est(key,self.cond_ind)
                #出力ファイルをest dictに登録する
                destf = self.make_tname('est',self.cond_est,'est')
                self.est_dict[key][self.index] = destf
                #estimate
                self.do_est(key,dat_file,destf,self.cond_est)
            self.cproc.exe()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Simulation for condition and plot roc and AUC')
    parser.add_argument('--tflag', '-t',default=False, action = "store_true",help='test flag not qsub')
    parser.add_argument('--condition', '-c',default='script/cond/emwer_test_condition.txt',type=str,help='condition filename')
    parser.add_argument('--hap_file', '-hp',default='tmp/test/hap_init_500.mimhap',type=str,help='condition filename')
    args = parser.parse_args()
    cp = cproc(False,args.tflag)
    cond_file = args.condition
    root_dir = "tmp/" + (cond_file.split(".")[0]).split("/")[-1] + "/"
    cdt = Conducter_sim(root_dir,cp,unit = 5)
    cdt.conducter(cond_file,args.hap_file,args.recomb_file)
    print(cproc(False))
