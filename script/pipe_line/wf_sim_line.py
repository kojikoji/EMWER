#-*- coding: utf-8 -*-
import random
class Sim_line:
    # do wf_sim
    def wf_sim_run(self,ofn):
        # make name of sync
        self.afs_file = self.dir_dict['afs'] + self.cond + '.afs'
        # make command and add parameter option for mimicree
        cmd = './build/default/WF_sim'
        cmd = ' '.join([cmd,'-n',self.state['allnum']])
        cmd = ' '.join([cmd,'-r',self.state['rep']])
        cmd = ' '.join([cmd,'-p',str(2*int(self.state['pop']))])
        cmd = ' '.join([cmd,'-s',self.state['slc']])
        cmd = ' '.join([cmd,'-f',self.state['dep']])
        cmd = ' '.join([cmd,'-g',self.state['gen']])
        cmd = ' '.join([cmd,'-d',self.state['disc']])
        if 'flc' in self.state:
            cmd = ' '.join([cmd,'-u',self.state['flc']])            
        if 'bd' in self.state:
            cmd = ' '.join([cmd,'-b',self.state['bd']])            
        cmd = ' '.join([cmd,'-o',ofn])
        cmd = ' '.join([cmd,'-a',self.afs_file])
        # cmd add
        self.cp.add(cmd)
   #do simulation
    def do_sim(self,ofn):
        self.wf_sim_run(ofn)
