class Condtion_controler:
#register dictionary of condition list(key is condition, and value is list of value
    #parse condition file and register it to condition dictionary
    def parse_cond_file(self,cond_file):
        f = open(cond_file)
        self.cond_dict = {}
        for line in f:
            line = line.strip()
            #if there is empty row, process not done 
            if(len(line)>3):
                #structure is key:val \t val \t val..
                key = line.split(':')[0]
                vals = line.split(':')[1].split('\t')
                #if top of row is # not registered
                if key[0] != '#':
                    self.cond_dict[key] = vals
    # make list of condtion in key-val_key-val_key-val
    def make_cond_list(self,cond_file='',setflag=False):
        if(cond_file != ''):
            self.parse_cond_file(cond_file)
        cond_list = ['']
        for key , vals in self.cond_dict.items():
            #make list of cond added for key's val
            new_cond_list = []
            for val in vals:
                elt = key+'-'+val
                for cond in cond_list:
                    #at first new condtion is key-val
                    if cond == '':
                        new_cond = elt
                    else:
                        #next key-val and key-val to key-val_key-val
                        new_cond = '_'.join([cond,elt])
                    new_cond_list.append(new_cond)
            cond_list = new_cond_list
        self.cond_list = cond_list
        if setflag:
            self.set_cond(cond_list[0])
    # change state (replicate selection) for cond
    def set_cond(self,cond):
        self.state = {}
        #condition is key-val_key-val_key-val
        self.cond = cond
        #[key-val,..]
        elts = cond.split('_')
        for elt in elts:
            #list [key,val]
            keyval = elt.split('-')
            self.state[keyval[0]] = keyval[1]
        if (not('genl' in self.state or 'repl' in self.state)):
            self.gen_rep_expander()
    # make genl and repl from rep and gen
    def gen_rep_expander(self):
        if not('gen' in self.state and 'rep' in self.state):
            print("#ERROR!!:condtion file is strange which contain only either gen or rep")
            exit(-1)
        uni_repl = map(str,range(0,int(self.state['rep'])))
        uni_genl = self.state['gen'].split(',')
        repl = []
        genl = []
        for rep in uni_repl:
            for gen in uni_genl:
                repl.append(rep)
                genl.append(gen)
        self.state['genl'] = ",".join(genl)
        self.state['repl'] = ",".join(repl)
