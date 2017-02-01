# -*- coding: utf-8 -*-
from cproc import cproc
import sys
def wfenc_estimate(self):
    #variables
    ifn = self.DDIR + "/sim" + "_s" + self.slc + ".dat "
    ifnn = self.DDIR + "/sim" + "_s" + self.slc + ".dat"
    ifn0 = self.DDIR + "/sim" + "_s0.dat "
    ln = int(self.Dnum)*2
    uni = 50
    init = 1
    ct = int(ln/uni) + 1
    for i in range(1,ct):
        end = init + uni-1
        if(end > ln):
            end = ln
        sfn = ifn + str(i)
        sfn = self.DDIR + "/tmp/sfnnc" + str(i) + ".dat "
        sfn0 = self.DDIR + "/tmp/sfnnc0" + str(i) + ".dat "
        srlt = self.DDIR + "/tmp/sncrlt" + str(i) + ".dat "
        cmdsed = "sed -n " + ",".join([str(init),str(end)]) + "p " + ifn + "> " + sfn
        cmdsed0 = "sed -n " + ",".join([str(init),str(end)]) + "p " + ifn0 + "> " + sfn0
        cmdest = "./bin/WF_test" \
               + " -p " + self.pop \
               + " -g " + self.gen \
               + " -s " + self.slc \
               + " -i " + sfn \
               + " -o " + srlt \
               + " -f " 
        self.cp.add(cmdsed)
        self.cp.add(cmdsed0)
        self.cp.add(cmdest)
        cmdrmsfn = "rm " + sfn + sfn0
        self.cp.add(cmdrmsfn)
        init += uni
        self.cp.mkqsf("ncestimate")
    self.cp.qsub("ncestimate","labeling")
    #conect file
    rlt = self.DDIR + "/rlt_wfenc.dat "
    crlt = self.DDIR + "/tmp/rlt_wfenc.dat "
    cmdrm = "rm -f " + rlt
    cmdtch = "touch " + rlt
    self.cp.add(cmdrm)
    self.cp.add(cmdtch)
    for i in range(1,ct):
        srlt = self.DDIR + "/tmp/sncrlt" + str(i) + ".dat "        
        cmdcat = "cat " + rlt + srlt + "> " + crlt
        cmdcp = "cp " + crlt + rlt
        self.cp.add(cmdcat)
        self.cp.add(cmdcp)
        cmdrmsrlt = "rm " + srlt
        self.cp.add(cmdrmsrlt)
    self.cp.mkqsf("ncconect")
    self.cp.qsub("ncconect","ncestimate")
    #java -Xmx1g -jar /home/ykojima/Projects/wfeed/script/mimicree/toysoft/MimicrEESummary.jar --output-format sync --threads 2 --recombination-rate /home/ykojima/Projects/wfeed/tmp/mimicree/dmel.rr.txt --output-mode 45 --additive /home/ykojima/Projects/wfeed/tmp/mimicree/sim_p1000_s0.01_r6/snps/s0.01_r6.txt --output-file /home/ykojima/yyProjects/wfeed/tmp/mimicree/sim_p1000_s0.01_r6/tmp/sim_g45.sync
    
