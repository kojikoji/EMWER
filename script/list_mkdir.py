#-*- coding: utf-8 -*-
import os
# if there is no directory, make it 
def safe_mkdir(dirname,rflag = False):
    if os.path.exists(dirname)!=True:
        os.mkdir(dirname)

class List_mkdir:
    def list_mkdir(self,root_dir='',tag='',list_dir_key=''):
        if root_dir != '':
            self.root_dir = root_dir
        # if tag is specified, add the tag to root dir name 
        if tag != '':
            self.root_dir = self.root_dir + '/' + tag + '/'
        if list_dir_key != '':
            self.list_dir_key = list_dir_key
        print(self.root_dir)
        safe_mkdir(self.root_dir,True)
        #initialize dictionary of directory
        self.dir_dict = {}
        self.tdir_dict = {}
        for key in self.list_dir_key:
            dir_name = self.root_dir + key + '/'
            self.dir_dict[key] = dir_name
            safe_mkdir(dir_name)
            #make tmp folder in each dir
            tdir_name = self.root_dir + key + '/tmp/'
            self.tdir_dict[key] = tdir_name
            safe_mkdir(tdir_name)
