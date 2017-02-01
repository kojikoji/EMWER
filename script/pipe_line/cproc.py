# -*- coding: utf-8 -*-
import subprocess
import os.path
import sys
class cproc:
    def __init__(self,cflag=False,tflag=False):
        # if cflag is True, program counter is read from fname
        #if not exist, program counter is set 0
        self.fname = ".program_count_test"
        self.cmdlist = []
        self.pc=0
        self.jsize="2"
        self.jtype=""
        if os.path.isfile(self.fname) and cflag:
            pcin = open(self.fname,'r')
            self.pc = int(pcin.read())
            pcin.close()
        # number of undo command is initialized by program counter
        self.nud = self.pc
        self.qsc = dict()
        # if tflag is true, qsub wont be done
        self.tflag = tflag
    #adding command to self.cmdlist
    def add(self,cmd,*depl):
        if len(depl) > 0:
            cmd = cmd + " " + " ".join(depl)
        if self.nud > 0:
            #if number of undo command is remain, dont register cmd 
            self.nud -= 1
            #print(pc)
        else:
            #if nud is 0, register cmd
            self.cmdlist.append(cmd)
    #execute command in self.cmdlist
    def exe(self,pflag=True):
        flag = 0
        for cmd in self.cmdlist:
            if pflag :
                print(cmd)
            flag = 0
            if self.tflag == False:
                flag =subprocess.call(cmd,shell=True)
            if flag != 0:
                #if error is occured, print out the command 
                print("#Error @ exe!!: %s" % cmd)
                break
            else:
                # increase prgram counter 
                self.pc += 1
        pcout = open(self.fname,'w')
        pcout.write(str(self.pc))
        pcout.close()
        self.cmdlist = []
        if flag != 0:
            exit(1)
    def fnget(self,fulfn):
        lst = fulfn.split("/")
        ans =lst[len(lst)-1]
        return(ans)
    def mkqsf(self,name="qsf"):
        name = self.fnget(name)
        if os.path.exists(".qsub")!=True:
            os.mkdir(".qsub")
        if os.path.exists(".qlog")!=True:
            os.mkdir(".qlog")
        if name not in self.qsc:
            self.qsc[name] = 0
        #qsc th file of command for qsub
        self.qsc[name] += 1
        qsfn = ".qsub/" + name + "." + str(self.qsc[name])
        qsf = open(qsfn,'w')
        #cmdlist to qsfstr
        qsfstr = "\n".join(self.cmdlist)
        qsf.write(qsfstr)
        qsf.close()
        self.cmdlist = []
    def job_size(self,size):
        self.jsize = str(size)
    def job_type(self,tp):
        if tp != "mjob":
            self.jtype = tp + ","
    def qsub(self,name="qsf",*depl):
        #first delete pre logs
        cmd = "rm .qlog/"+name+"*"
        subprocess.call(cmd,shell=True)
        name = self.fnget(name)
        tdepl = []
        for dep in depl:
            tdepl.append(self.fnget(dep))
        depl = tdepl
        #make the content of qsub file
        qsarstr = "#!/bin/bash\n" \
                  + "#$ -S /bin/bash\n" \
                  + "#$ -cwd\n" \
                  + "#$ -M tiisaishima@gmail.com\n" \
                  + "#$ -l "+self.jtype+"s_vmem="+self.jsize+"G,mem_req="+self.jsize+"G\n" \
                                                 + "bash .qsub/" + name + ".$SGE_TASK_ID "
        #make qsub file qsar.sh 
        qsarf = open(".qsub/"+name+"ar.sh",'w')
        qsarf.write(qsarstr)
        qsarf.close()
        anum_str = "1-" + str(self.qsc[name])
        #reset
        self.qsc[name] = 0
        self.jsize="2"
        self.jtype=""
        cmd = "qsub " + " -t " + anum_str + " -N " + name
        cmd +=  " -o ./.qlog/"
        cmd += " -e ./.qlog/"
        #add dependency
        if len(depl) > 0:
            cmd += " -hold_jid "
            cmd += ",".join(depl)
        cmd += " .qsub/"+name+"ar.sh"
        #
        #array job is released
        flag = 0
        if self.tflag==False :
            flag = subprocess.call(cmd,shell=True)
        else:
            print(cmd)
        if flag != 0:
            #if error is occured, print out the command 
            print("#Error @ qsub!! : %s" % cmd)
            exit(1)

        
if __name__=="__main__":
    #continuous flag is true when command line is c
    cflag = sys.argv[len(sys.argv)-1] == "c"
    cp = cproc(tflag = True)
    # pc mean number of command undo
    # init_pc mean reach point
    cp.add("ls -a")
    cp.add("ls ..")
    cp.add("ls")
    cp.exe()
    cp.add("ls",'-a','-l')
    cp.exe()
    #cp.add_args(['-a','-l'])    
    #cp.mkqsf()
    cp.add("ls -a")
    cp.add("less test.txt")
    cp.add("ls -a")
    #cp.job_type("sjob")
    cp.mkqsf()
    cp.qsub()    
    
