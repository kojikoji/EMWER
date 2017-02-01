# -*- coding: utf-8 -*-
import random as rd
import sys
def make_oneloci(num,thresh):
    ans = ""
    ac =0
    cc =0
    for i in range(num):
        val = rd.random()
        if i%2==0 and i!=0:
            ans += " "
        if val > thresh:
            ans += "A"
            ac+=1
        else:
            ans += "C"
            cc+=1
    if ac*cc == 0:
        ans =""
    return ans

if __name__ == "__main__":
    lnum = int(sys.argv[1])
    num = int(sys.argv[2])
    if len(sys.argv) > 3 :
        filename = sys.argv[3]
    else:
        filename = "hap"
        filename += str(num)
        filename += ".mimhap"

    data = ""
    for i in range(lnum):
        line = ""
        while line == "":
            thresh = rd.random()
            line = make_oneloci(num,thresh)
        data += line
        data += "\n"
    out_f = open(filename,'w')
    out_f.write(data)
    out_f.close()
