#-*- coding: utf-8 -*-
import random
class Sim_line:
    def sample_hap(self):
        #extract part of haplotype and to list
        #make hap file
        self.hap_file = self.dir_dict['hap'] + self.cond + '.mimhap'
        cmd = "python script/sample_hap.py"
        cmd = " ".join([cmd,"-p",self.state['pop']])
        cmd = " ".join([cmd,"-i",self.ohap_file])
        cmd = " ".join([cmd,"-o",self.hap_file])
        self.cp.add(cmd)
    #decide selected snps
    def make_selc_snp(self):
        # make name of selected snp file
        self.slc_file = self.dir_dict['snp'] + self.cond + '.snps'
        # make command and add parameter option for add snp
        cmd = 'python script/mimicree/toysoft/pick-random-addtive-snps.py '
        cmd = ' '.join([cmd,'-n',self.state['slcnum']]) ## !!!!!
                       #dominance
        cmd = ' '.join([cmd,'-e','0.5'])
        cmd = ' '.join([cmd,'-m','0.9'])
        cmd = ' '.join([cmd,'-s',self.state['slc']])
        cmd = ' '.join([cmd,'--loci-count',str(self.loc_num)])
        cmd = ' '.join([cmd,'--input',self.hap_file])
        cmd = ' '.join([cmd,'>',self.slc_file])
        self.cp.add(cmd)
    # do mimicree
    def mimicree_run(self):
        # make name of sync
        self.presync_file = self.tdir_dict['sim'] + self.cond + '.sync'
        # make command and add parameter option for mimicree
        cmd = 'java -Xmx1g -jar  script/mimicree/toysoft/MimicrEESummary.jar --output-format sync --threads 1'
        cmd = ' '.join([cmd,'--haplotypes-g0',self.hap_file])
        cmd = ' '.join([cmd,'--recombination-rate',self.recomb_file])
        # generations of output  
        cmd = ' '.join([cmd,'--output-mode',self.state['gen']])
        cmd = ' '.join([cmd,'--replicate-runs',self.state['rep']])
        cmd = ' '.join([cmd,'--additive',self.slc_file])
        cmd = ' '.join([cmd,'--output-file',self.presync_file])
        # cmd add
        self.cp.add(cmd)
    #poisson sample
    def poisson_sample(self):
        self.sync_file = self.dir_dict['sim'] + self.cond + '.sync'
        cmd = 'python script/mimicree/poisson-3fold-sample.py '
        cmd = ' '.join([cmd,'--input',self.presync_file]) ## !!!!!
        cmd = ' '.join([cmd,'--coverage',self.state['dep']]) ## !!!!!
        cmd = ' '.join([cmd,'>',self.sync_file])
        self.cp.add(cmd)
    #make labels for selected snp
    def make_label(self):
        # make name of label
        self.label_file = self.dir_dict['label'] + self.cond + '.lb'
        # make command and add parameter option for label
        cmd = 'python script/mimicree/toysoft/rocr-generate-labellist.py'
        cmd = ' '.join([cmd,self.sync_file])
        cmd = ' '.join([cmd,self.slc_file])
        cmd = ' '.join([cmd,'>',self.label_file])    
        self.cp.add(cmd)
   #do simulation
    def do_sim(self):
        #self.sample_hap()
        self.make_selc_snp()
        self.mimicree_run()
        self.poisson_sample()
        self.make_label()
        self.cp.job_size("8")
