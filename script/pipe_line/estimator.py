#-*- coding: utf-8 -*-
class Estimator:
    # est file を keyに応じて格納するdictを作る
    def est_init(self,keys):
       self.est_dict = dict()
       for key in keys:
           #ファイル名を推定回数分各keyについて用意
           est_list = [""]*self.est_iter
           self.est_dict[key] = est_list
    def set_est(self,key,tag):
        self.cond_est = tag + '_est-'+key
    def do_est(self,key,ifn,ofn,tag):        
        #キーに応じてコマンドが変わる
        cmd = self.wfe_command(ifn,ofn)
        if key == 'wfe':
            cmd = self.wfe_command(ifn,ofn)
        #キーに応じてコマンドが変わる
        if key == 'cmha':
            cmd = self.cmh_command(ifn,ofn)
        #キーに応じてコマンドが変わる
        if key == 'cmh':
            cmd = self.cmhtp_command(ifn,ofn)
        #キーに応じてコマンドが変わる
        if key == 'gcmh':
            cmd = self.gcmh_command(ifn,ofn)
        #キーに応じてコマンドが変わる
        if key == 'mlwf':
            cmd = self.mlwf_command(ifn,ofn)
        #キーに応じてコマンドが変わる
        if key == 'bbgp':
            cmd = self.bbgp_command(ifn,ofn)
        self.cp.add(cmd)
    def concat_est(self,key,ofn):
        cmd = "rm " + ofn
        self.cp.add(cmd)
        for dest in self.est_dict[key]:
            cmd = "awk \'{print $0}\'"
            cmd = " ".join([cmd ,dest])
            cmd = " ".join([cmd,">>",ofn])
            self.cp.add(cmd)

        #self.cp.add(cmdp)
        #self.cp.add(cmd)
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
    #wfeのコマンド作成
    def wfe_command(self,ifn,ofn):
        cmd = "./build/default/WF_main "
        cmd = " ".join([cmd,"-d",self.state['disc']])
        cmd = " ".join([cmd,"-p",self.state['pop']])
        cmd = " ".join([cmd,"-i",ifn])
        cmd = " ".join([cmd,"-o",ofn])
        # f flag optimize selection
        if self.state['opt'] == "slc":
            cmd = " ".join([cmd,"-f"])
        # l flag optimize population
        else:
            cmd = " ".join([cmd,"-l"])
        # one estimate for all data
        if 'all' in self.state:
            cmd = " ".join([cmd,"-a"])            
        if 'afs' in self.state:
            if self.state['afs'] == 't':
                cmd = " ".join([cmd,"-q",self.afs_file])
        return(cmd)
    #wfeのコマンド作成
    def mlwf_command(self,ifn,ofn):
        cmd = "python script/ml_wf/threelocus/wf_ml.py"
        cmd = " ".join([cmd,"-l",self.state['hapnum']])
        cmd = " ".join([cmd,"-g",self.state['gen']])
        cmd = " ".join([cmd,"-p",self.state['pop']])
        cmd = " ".join([cmd,"-i",ifn])
        cmd = " ".join([cmd,"-o",ofn])
        # one estimate for all data
        if 'afs' in self.state:
            if self.state['afs'] == 't':
                cmd = " ".join([cmd,"-a",self.afs_file])
        return(cmd)
    #cmhのコマンド作成
    def cmhtp_command(self,ifn,ofn):
        cmd = "Rscript"
        cmd = " ".join([cmd,"script/mimicree/run_cmhtp.r"])
        cmd = " ".join([cmd,ifn])
        cmd = " ".join([cmd,ofn])
        return(cmd)
    #cmhのコマンド作成
    def cmh_command(self,ifn,ofn):
        cmd = "Rscript"
        cmd = " ".join([cmd,"script/mimicree/run_cmh_all.r"])
        cmd = " ".join([cmd,ifn])
        cmd = " ".join([cmd,ofn])
        return(cmd)
    #cmhのコマンド作成
    def gcmh_command(self,ifn,ofn):
        cmd = "Rscript"
        cmd = " ".join([cmd,"script/mimicree/run_gcmh.r"])
        cmd = " ".join([cmd,ifn])
        cmd = " ".join([cmd,ofn])
        return(cmd)
    def bbgp_command(self,ifn,ofn):
        cmd = "Rscript"
        cmd = " ".join([cmd,"script/BBGP/bbgp_do.r"])
        cmd = " ".join([cmd,ifn])
        cmd = " ".join([cmd,ofn])
        return(cmd)
