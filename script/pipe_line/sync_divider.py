import os.path
class Sync_divider:            
    def get_tag(self,fname):
        return(fname.split('/')[-1].split('.')[0])
    # change the index of estimate is set for i
    def divide_file(self,ifn,ofn):
        cmd = 'python script/index_divider.py'
        cmd = " ".join([cmd,'-i',ifn])
        cmd = " ".join([cmd,'-o',ofn])
        cmd = " ".join([cmd,'-x',str(self.index)])
        cmd = " ".join([cmd,'-d',str(self.est_iter)])            
        self.cp.add(cmd)
     #devide file for index th estimate
    def divide_file_dat(self,ifn,ofn):
        self.condind = self.cond + '_n-' + str(self.index)
        self.dat_file = self.tdir_dict['dat'] + self.condind + '.dat'
        unit = self.unit
        # sed -n st,edp sync dsync
        cmd = 'python script/index_divider.py'
        cmd = " ".join(['-i',self.all_dat_file])
        cmd = " ".join(['-o',self.dat_file])
        cmd = " ".join(['-x',self.index])
        cmd = " ".join(['-d',self.est_iter])            
        self.cp.add(cmd)
