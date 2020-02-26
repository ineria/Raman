# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import matplotlib.pyplot as plt
import numpy_indexed as npi
from math import floor

"""This code takes 2D map ratio data and transforms them into a 1D plot"""

data_dir="" #specify data directory

"""
    In this part of the code the parameters are defined. This specifies where the data are found
    Takes filenames and names for components
    Pay attention to enter correct STEPSIZE !!!
"""
date= "date/"
filename="name"
filename_small="name_component1_wavenumber.txt"
filename_large="name_component2_wavenumber.txt"
large= "component2"
small= "component1"

step=20 ##stepsize in micrometer between each spectrum

"""Functions to be used in data analysis later"""
def average_with_window(data, window = 5):
  """
    data is a flat array
    window is the number of entries from the array that shall be averaged
    returns an array of length floor(len(data)/window) with average values.
  """
  avgs = []
  for i in range(floor(len(data)/window)):
    start_of_window = i * 5
    avgs.append(np.average(data[start_of_window : start_of_window + window]))
    
  return avgs


def std_with_window(data, window = 5):

  std = []
  for i in range(floor(len(data)/window)):
    start_of_window = i * 5
    std.append(np.std(data[start_of_window : start_of_window + window]))
    
  return std


data_small=np.loadtxt(data_dir+date+filename_small, dtype='float', comments='#', delimiter='\t').T                   
x=data_small[0]
I_small=data_small[2]

data_large=np.loadtxt(data_dir+date+filename_large, dtype='float', comments='#', delimiter='\t').T                    
I_large=data_large[2]

x, I_small_mean=npi.group_by(x).mean(I_small)
x, I_large_mean=npi.group_by(x).mean(I_large)
depth=np.asarray(np.arange(0,len(I_small_mean))*step*np.sin(2*np.pi/180))

"""1D plot of ratio value over position averaged in columns"""
plt.plot(depth[::-1],(I_small_mean/(I_large_mean+I_small_mean)),'k+')
plt.ylabel(r'$\frac {I_{'+ small +'}}{I_{'+ large +'}+I_{'+ small +'}}}$',fontsize=16)
plt.xlabel(r'$depth\ below\ surface\ [\mu m]$')
#plt.ylim(0.2,0.8)
plt.savefig(data_dir+"RatiomapAnalysis/"+filename+'1D_fraction_plot_zoom.svg', transparent=True)
plt.show()


"""average of average to reduce data (average over 5 data points)"""

plt.errorbar(average_with_window(depth, 5)[::-1],average_with_window((I_small_mean/(I_large_mean+I_small_mean)), 5), std_with_window((I_small_mean/(I_large_mean+I_small_mean)), 5), linestyle='None', marker='o')
plt.ylabel(r'$\frac {I_{'+ small +'}}{I_{'+ large +'}+I_{'+ small +'}}}$',fontsize=16)
plt.xlabel(r'$avg.\ depth\ below\ surface\ [\mu m]$')
#plt.ylim(0.2,0.8)
plt.savefig(data_dir+"RatiomapAnalysis/"+filename+'1D_fraction_plot_avg_zoom.svg', transparent=True)
plt.show()

print(len(average_with_window((I_small_mean/(I_large_mean+I_small_mean)), 5)))
print(np.average(average_with_window((I_small_mean/(I_large_mean+I_small_mean)), 5)))
