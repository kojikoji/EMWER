#-*- coding: utf-8 -*-
class Hogehoge:
    def __init__(self):
        self.value = 0
    def plus(self,value):
        return(self.value + value)
import random
class Sim_line:
    #decide selected snps
    def make_selc_snp(self,ifn,ofn):
        # make name of selected snp file
        self.slc_file = self.dir_dict['snp'] + self.cond + '.snps'
        # make command and add parameter option for add snp
        cmd = 'python script/mimicree/toysoft/pick-random-addtive-snps.py '
        cmd = ' '.join([cmd,'-n',self.state['slcnum']]) ## !!!!!
                       #dominance
        cmd = ' '.join([cmd,'-e','0.5'])
        cmd = ' '.join([cmd,'-m','0.9'])
        cmd = ' '.join([cmd,'-s',self.state['slc']])
        cmd = ' '.join([cmd,'--loci-count',self.state['num']])
        cmd = ' '.join([cmd,'--input',ifn])
        cmd = ' '.join([cmd,'>',ofn])
        self.cp.add(cmd)
    # do mimicree
    def mimicree_run(self,ifn,ofn,snpfn):
        # make name of sync
        self.presync_file = self.tdir_dict['sim'] + self.cond + '.sync'
        # make command and add parameter option for mimicree
        cmd = 'java -Xmx1g -jar  script/mimicree/toysoft/MimicrEESummary.jar --output-format sync --threads 1'
        cmd = ' '.join([cmd,'--haplotypes-g0',ifn])
        cmd = ' '.join([cmd,'--recombination-rate',self.recomb_file])
        # generations of output  
        cmd = ' '.join([cmd,'--output-mode',self.state['gen']])
        cmd = ' '.join([cmd,'--replicate-runs',self.state['rep']])
        cmd = ' '.join([cmd,'--additive',snpfn])
        cmd = ' '.join([cmd,'--output-file',ofn])
        # cmd add
        self.cp.add(cmd)
    #poisson sample
    def poisson_sample(self,ifn,ofn):
        cmd = 'python script/mimicree/poisson-3fold-sample.py '
        cmd = ' '.join([cmd,'--input',ifn]) ## !!!!!
        cmd = ' '.join([cmd,'--coverage',self.state['dep']]) ## !!!!!
        cmd = ' '.join([cmd,'>',ofn])
        self.cp.add(cmd)
    #make labels for selected snp
    def make_label(self,ifn,ofn,snpfn):
        # make name of label
        # make command and add parameter option for label
        cmd = 'python script/mimicree/toysoft/rocr-generate-labellist.py'
        cmd = ' '.join([cmd,ifn])
        cmd = ' '.join([cmd,snpfn])
        cmd = ' '.join([cmd,'>',ofn])    
        self.cp.add(cmd)
   #do simulation
    def do_sim(self,ifn,ofn,tag):
        lbfn=self.make_name('label',tag,'lb')
        #self.sample_hap()
        snpfn = self.make_name('snp',tag,'snp')
        self.make_selc_snp(ifn,snpfn)
        mimfn = self.make_tname('sim',tag,'psync')
        self.mimicree_run(ifn,mimfn,snpfn)
        self.poisson_sample(mimfn,ofn)
        self.make_label(mimfn,lbfn,snpfn)
        self.cp.job_size("8")
