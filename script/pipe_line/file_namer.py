# -*- coding: utf-8 -*-
class File_namer:
    def make_name(self,dirtag,tag,ext):
        name = self.dir_dict[dirtag] + tag + "." + ext
        return(name)
    def make_tname(self,dirtag,tag,ext):
        name = self.tdir_dict[dirtag] + tag + "." + ext
        return(name)
