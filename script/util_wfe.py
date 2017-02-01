#!/usr/env python
#-*- coding: utf-8 -*-
from numpy.random import *
import numpy as np
import sys
import random
import matplotlib.pyplot as plt
def sum_regulalize(numpy_array,regulalized_sumed_value):
    sumed_value = sum(numpy_array)
    print(sumed_value)
    regulalized_numpy_array = (regulalized_sumed_value/sumed_value)*numpy_array
    return(regulalized_numpy_array)
def cumrative(numarray):
    cumrative_array = np.zeros(len(numarray))
    cumrative_array[0] = numarray[0]
    for i in range(1,len(numarray)):
        cumrative_array[i] = cumrative_array[i-1] + numarray[i]
    return(cumrative_array)
def cumrative_sample(cumrative_array,pointer):
    index = 0
    for cumrative_value in cumrative_array:
        if cumrative_value > pointer:
            break;
        index += 1
    return(index)
class Bin:
    def __init__(self,low,high,scale):
        self.low = float(low)
        self.high = float(high)
        self.scale = float(scale)
    def get(self,value):
        scaled_value = ((self.high - self.low)/self.scale)*value
        if(scaled_value < 0 or scaled_value > (self.high - self.low)):
            print >> sys.stderr,'Error:This value is bigger than upper bound'
        value_in_bin = self.low + scaled_value
        return value_in_bin
class Spectrum_uniform:
    def set_low_high(self,low,high):
        self.low = low
        self.high = high
    def set_spectrum(self,spectrum):
        lower_bound = self.low
        self.bin_width = (self.high - self.low)/float(len(spectrum))
        self.cumrative_spectrum = cumrative(spectrum)
        self.cumrative_spectrum_anterior = self.cumrative_spectrum - spectrum
        print(self.cumrative_spectrum)
        self.bin_list = list()
        for measure in spectrum:
            upper_bound = lower_bound + self.bin_width
            self.bin_list.append(Bin(lower_bound,upper_bound,measure))
            lower_bound = upper_bound
    def sample(self,low,high):
        pointer = random.uniform(0,self.cumrative_spectrum[-1])
        sample_index = cumrative_sample(self.cumrative_spectrum,pointer)
        inner_pointer = pointer - self.cumrative_spectrum_anterior[sample_index]
        sample_value = self.bin_list[sample_index].get(inner_pointer)
        if(sample_value < low or sample_value > high):
            sample_value = self.sample(low,high)
        return(sample_value)
if __name__ == "__main__":
    su = Spectrum_uniform()
    su.set_low_high(0,1)
    data = np.loadtxt("tmp/afs_default.afs")
    su.set_spectrum(data)
    sulist = list()
    for i in range(10000):
        sulist.append(su.sample(0.1,0.9))
    plt.figure(figsize=(8,6),dpi=80)
    plt.subplot(1,1,1)
    plt.hist(sulist,bins=100)
    plt.show()
