# -*- coding: utf-8 -*-


import os
import matplotlib.pyplot as plt
import numpy as np
import glob
from matplotlib.ticker import (MultipleLocator)
from scipy import stats

## setup files etc
#data 
data_dir = "" #specify data directory
data_file_count = int(len(os.listdir(data_dir))+1) # +1 later!
data_files =sorted(glob.glob(data_dir+'*Average*_Acq*.txt'))

#set compounds and filename
name="filename"
large= "component 1"
small= "component 2"


if 'small=component 2':
        wave_small,wave_large=(247,730) #lines in txt file where relevant peak data can be found
        
elif 'small=component 3':
        wave_small,wave_large=(622,189)

## selecting specific wavenumber from files and extract intensity for each spectrum to new file

def get_component_depth_profile (w1,w2):##w1,w2 wavenumber of interest for two compounds
    for i,filename in enumerate(data_files):
        np.loadtxt(filename, delimiter='\t',comments='#')   
        
        with open(data_dir+'temp1.txt', 'w') as outfile:
            for filename in data_files:
                data_list = open(filename).readlines()
                outfile.write(data_list[w1])

        with open(data_dir+'temp2.txt', 'w') as outfile:
            for filename in data_files:
                data_list = open(filename).readlines()
                outfile.write(data_list[w2])

    with open(data_dir+'peakdata.txt', 'w') as file3:
        with open(data_dir+'temp1.txt', 'r') as file1:
            with open(data_dir+'temp2.txt', 'r') as file2:
                for line1, line2 in zip(file1, file2):
                    print(line1.strip(), line2.strip(), file=file3) 
                    
    os.remove(data_dir+'temp1.txt')
    os.remove(data_dir+'temp2.txt')



              
get_component_depth_profile(wave_small,wave_large) 
             
## plot profile

int_small,int_large=np.loadtxt(data_dir+'peakdata.txt').T[[1,3]]
x=[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0]
#x=[0,5,10,15,20,25,30,33,50,75,100]
fig, ax1 =plt.subplots()

ax1.plot(x,int_small/(int_small+int_large), 'o', color='#003e7e',markersize=10)
ax1.set_ylabel(r'$\frac {I_{'+ small +'}}{I_{'+ large +'}+I_{'+ small +'}}$', color='#003e7e', fontsize=16)
#ax1.set_ylim(0,1)
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
ax1.set_xlabel(r'$amount\ of\ $' + small + r'$[\%]$', color='#003e7e', fontsize=16)

slope, intercept, r_value, p_value, std_err = stats.linregress(x,int_small/(int_small+int_large))               
plt.plot(x,intercept+slope*np.array(x),'r--',label='y={:.4f}x+{:.4f}'.format(slope,intercept))


plt.legend(fontsize=16)               
plt.grid(True, which='both')
plt.savefig(data_dir+name+'_calibration_v2.svg', Transparent=True)
os.remove(data_dir+'peakdata.txt')

