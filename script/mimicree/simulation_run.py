# -*- coding: utf-8 -*-
from cproc.cproc import cproc
cp = cproc(False)
import sys
from wfeed_estimate import *
from wfenc_estimate import *
from bbgp_estimate import *
from cmh_estimate import *
class sim_run:
    def __init__(self):
        self.slc = "0.02"
        #variables
        self.WDIR="/home/ykojima/Projects/wfeed/tmp/mimicree"
        self.SDIR="/home/ykojima/Projects/wfeed/script/mimicree"
        self.gen = "15,30,45,60"
        self.pop = "1000"
        self.fld = "100"
        self.Dnum = "1000"
        self.dnum = "12"
        self.DDIR = self.WDIR + "/sim" \
                    + "_p" + self.pop \
                    + "_gen" + self.gen \
                    + "_slc" + self.slc \
                    + "_f" + self.fld \
                    + "_dnum" + self.dnum \
                    + "_Dnum" + self.Dnum
        self.cp = cp
    def do_setting(self):
    #command for setting
        cmdmkdir = "mkdir " + self.DDIR + " " + self.DDIR + "/tmp " + self.DDIR + "/snps" 
        self.cp.add(cmdmkdir)
        self.cp.mkqsf("setting")
        self.cp.qsub("setting")
    def make_simfile(self,slc):
        fnsim = self.DDIR + "/simi" + "_s" + slc + ".dat "
        fnsimc = self.DDIR + "/simc" + "_s" + slc + ".dat "
        fnrltb = "_p" + self.pop \
                 + "_gen" + self.gen \
                 + "_slc" + slc \
                 + "_f" + self.fld \
                 + "_dnum" + self.dnum \
                 + "_Dnum" + self.Dnum \
                 + ".dat"
        cmdsim = "./bin/WF_sim" \
                 + " -N " + self.Dnum \
                 + " -n " + self.dnum \
                 + " -p " + self.pop \
                 + " -g " + self.gen \
                 + " -f " + self.fld \
                 + " -s " + slc \
                 + " -o " + fnsim \
                 + " -c " + fnsimc
        self.cp.add(cmdsim)
        self.cp.mkqsf("simulation"+slc)
        self.cp.qsub("simulation"+slc,"setting")
    def do_sim(self):
        self.make_simfile(self.slc)
        self.make_simfile("0")
    def make_labels(self):
        TAG=self.DDIR.replace('/','_')
        self.WDIR="/home/ykojima/Projects/wfeed/tmp/mimicree"
        SDIR="/home/ykojima/Projects/wfeed/script/mimicree"
        fnsim = self.DDIR + "/simi" + "_s" + self.slc + ".dat "
        fnsimc = self.DDIR + "/simc" + "_s" + self.slc + ".dat "
        fnsim0 = self.DDIR + "/simi" + "_s" + "0" + ".dat "
        fnsimc0 = self.DDIR + "/simc" + "_s" + "0" + ".dat "
        clbfn = self.DDIR + "/tmp/lb_s" + self.slc + "_r6.txt "
        lbfn = self.DDIR + "/lb_s" + self.slc + "_r6.txt "
        ofnda = self.DDIR + "/sim" + "_s" + self.slc + ".dat "
        ofnda0 = self.DDIR + "/sim" + "_s" + "0" + ".dat "
        cmdawk1 = "awk \'{print 1}\' " + fnsim + "> " + clbfn
        cmdawk0 = "awk \'{print 0}\' " + fnsim0 + ">> " + clbfn
        cmdcat = "cat " + fnsim + "> " + ofnda
        cmdcat0 = "cat " + fnsim0 + ">> " + ofnda
        cmdcatc = "cat " + fnsimc + "> " + ofnda0
        cmdcatc0 = "cat " + fnsimc0 + ">> " + ofnda0
        self.cp.add(cmdawk1)
        self.cp.add(cmdawk0)
        self.cp.add(cmdcat)
        self.cp.add(cmdcat0)
        self.cp.add(cmdcatc)
        self.cp.add(cmdcatc0)
        self.cp.mkqsf("labeling")
        self.cp.qsub("labeling","simulation"+self.slc,"simulation"+"0")
sim_run.wfe_est = wfeed_estimate                
sim_run.wfenc_est = wfenc_estimate                
sim_run.bbgp_est = bbgp_estimate                
sim_run.cmh_est = cmh_estimate                
if __name__=="__main__":
    #variables
    #do simulation
    sr = sim_run()
    sr.do_setting()
    sr.do_sim()
    sr.make_labels()
    sr.wfe_est()
    sr.wfenc_est()
    sr.bbgp_est()
    sr.cmh_est()
    #java -Xmx1g -jar /home/ykojimaProjects/wfeed/script/mimicree/toysoft/MimicrEESummary.jar --output-format sync --threads 2 --recombination-rate /home/ykojimaProjects/wfeed/tmp/mimicree/dmel.rr.txt --output-mode 45 --additive /home/ykojimaProjects/wfeed/tmp/mimicree/sim_p1000_s0.01_r6/snps/s0.01_r6.txt --output-file /home/ykojimaProjects/wfeed/tmp/mimicree/sim_p1000_s0.01_r6/tmp/sim_g45.sync
    
