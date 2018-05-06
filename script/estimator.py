#-*- coding: utf-8 -*-
class Estimator:
    def est_init(self,keys):
       self.est_dict = dict()
       for key in keys:
           est_list = [""]*self.est_iter
           self.est_dict[key] = est_list
    def set_est(self,key,tag):
        self.cond_est = tag + '_est-'+key
    def do_est(self,key,ifn,ofn,tag):        
        cmd = self.wfe_command(ifn,ofn)
        if key == 'wfe':
            cmd = self.wfe_command(ifn,ofn)
        self.cp.add(cmd)
    def concat_est(self,key,ofn):
        cmd = "rm " + ofn
        self.cp.add(cmd)
        for dest in self.est_dict[key]:
            cmd = "awk \'{print $0}\'"
            cmd = " ".join([cmd ,dest])
            cmd = " ".join([cmd,">>",ofn])
            self.cp.add(cmd)
    def sync_to_dat(self,ifn,ofn,tag):
        #make dat file and afs file and ref_file
        self.ref_file = self.make_name('ref',tag,'ref')
        cmd = "python script/sync_to_dat_ref_base.py"
        cmd = " ".join([cmd,"-g",self.state['genl']])
        cmd = " ".join([cmd,"-r",self.state['repl']])
        cmd = " ".join([cmd,"-i",ifn])
        cmd = " ".join([cmd,"-o",ofn])
        cmd = " ".join([cmd,"-f",self.ref_file])
        if 'afs' in self.state:
            self.afs_file = self.make_name('afs',tag,'afs')
            cmd = " ".join([cmd,"-a",self.afs_file])
        if 'rbd' in self.state:
            cmd = " ".join([cmd,"-b",self.state['rbd']])
        self.cp.add(cmd)
    def wfe_command(self,ifn,ofn):
        cmd = "./build/default/WF_main "
        cmd = " ".join([cmd,"-d",self.state['disc']])
        cmd = " ".join([cmd,"-p",self.state['pop']])
        cmd = " ".join([cmd,"-i",ifn])
        cmd = " ".join([cmd,"-o",ofn])
        #which parameter estimate
        cmd = " ".join([cmd,"-f",self.state['opt']])        
        #which parameter estimate
        cmd = " ".join([cmd,"-t",self.state['thread']])        
        # one estimate for all data
        if 'all' in self.state:
            cmd = " ".join([cmd,"-a"])            
        # define selection for variant allele
        if self.cond_checker("allele","variant"):
            cmd = " ".join([cmd,"-v"])
        return(cmd)
