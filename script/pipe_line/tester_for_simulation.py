import unittest
import os
from conduct_simulation import Conducter_sim

class Conducter_sim_test(unittest.TestCase):
    def setUp(self):
        self.cdt = Conducter_sim('tmp/test/sim_root/')
        self.cdt.list_mkdir()
    def test_list_mkdir_there(self):
        self.cdt.list_mkdir()
        #there is dir for wfe
        self.assertEqual(True,os.path.exists('tmp/test/sim_root/wfe/'))
        self.assertEqual(True,os.path.exists('tmp/test/sim_root/cmh/tmp/'))
    def test_list_mkdir_dict(self):
        self.cdt.list_mkdir()
        #there is dir for wfe
        self.assertEqual('tmp/test/sim_root/hap/',self.cdt.dir_dict['hap'])
        self.assertEqual('tmp/test/sim_root/sim/tmp/',self.cdt.tdir_dict['sim'])
    def test_parse_cond_file(self):
        self.cdt.parse_cond_file('tmp/test/test_cond_sim.txt')
        self.assertEqual(['2'],self.cdt.cond_dict['rep'])        
        self.assertEqual(['0.05','0.1'],self.cdt.cond_dict['slc'])        
    def test_make_cond_list(self):
        self.cdt.make_cond_list('tmp/test/test_cond_sim.txt')
        self.assertEqual(2,len(self.cdt.cond_list))
    def test_set_cond(self):
        self.cdt.set_cond('a-3_box-188')
        self.assertEqual('188',self.cdt.state['box'])
        self.assertEqual('3',self.cdt.state['a'])
    def test_loc_number(self):
        self.cdt.hap_file = 'tmp/test/0-18.txt'
        self.cdt.sync_file = 'tmp/test/0-18.txt'
        self.cdt.divide_init()
        self.assertEqual(19,self.cdt.loc_num)
        self.assertEqual(7,self.cdt.est_iter)
        self.assertEqual(3,self.cdt.unit)
    def test_divide_file(self):
        self.cdt.list_mkdir()
        self.cdt.hap_file = 'tmp/test/0-18.txt'
        self.cdt.sync_file = 'tmp/test/0-18.txt'
        self.cdt.divide_init()
        self.cdt.set_index(self.cdt.est_iter-1)
        self.cdt.cond = '0-18'
        self.cdt.divide_file()
        self.assertEqual('tmp/test/sim_root/sim/tmp/0-18_n-6.sync',self.cdt.dsync_file)
        self.cdt.cp.exe()
        self.assertEqual('18\n',open('tmp/test/sim_root/sim/tmp/0-18_n-6.sync').read())
        self.cdt.set_index(3)
        self.cdt.divide_file()
        self.cdt.cp.exe()
        self.assertEqual('9\n10\n11\n',open('tmp/test/sim_root/sim/tmp/0-18_n3.sync').read())
    def test_do_sim(self):
        self.cdt.hap_file = 'tmp/test/hap_init_500.mimhap'
        self.cdt.recomb_file = 'tmp/test/dmel.rr.txt'
        self.cdt.make_cond_list('tmp/test/test_cond_sim.txt')
        self.cdt.list_mkdir()
        for cond in self.cdt.cond_list:
            print(cond)
            self.cdt.set_cond(cond)
            #self.cdt.do_sim()
        #self.cdt.cp.exe()
    def test_sample_hap(self):
        self.cdt.ohap_file = 'tmp/test/hap_sampler_test.txt'
        self.cdt.set_cond('pop-3') 
        self.cdt.sample_hap()
        self.cdt.cp.exe()
    def test_do_est(self):
        self.cdt.sync_file = 'tmp/test/sync_test.sync'
        self.cdt.set_cond('disc-50_pop-500_slc-0.05_opt-slc_gen-0,30_rep-2')
        self.cdt.unit=3
        self.cdt.divide_init()
        est_keys = ['wfe','cmh','bbgp']
        self.cdt.est_init(est_keys)
        for i in range(0,self.cdt.est_iter):
            self.cdt.set_index(i)
            self.cdt.divide_file()
            for key in est_keys:
                self.cdt.do_est(key)
        for key in est_keys:
            print(key)
            self.cdt.concat_est(key)
        self.cdt.cp.exe()
    
# do test
if __name__=="__main__":
    unittest.main()
