# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:57:15 2020

@author: ms01106

"""
##These are libraries that are imported in order to make the program work. Best just to leave these as they are now
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy import signal 
import glob
import os
import re
import matplotlib.cm
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib.text import Annotation






"""This is where you set up all your data.
I tried to write the program in a way that you only need to make changes here once and don't need to scroll through the entire script"""

data_dir = "//surrey.ac.uk/personal/HS129/ms01106/.System/Desktop/Katrin Raman Scripts/" ##useful to set up the data directory so that your script doesn't have to be in the same folder as your data
data_files =sorted(glob.glob(data_dir+'2020.01_gr bl 50k 600s large 0d rs 2200_Acq??.txt')) ##to identify all the data you want to read in they have to have some common element in the naming. Best works naming series the same and add Acq00 at the end. This then makes it easy to read in all data at the same time in the correct order
data_file_count = int(len(os.listdir(data_dir))+1)
series='2020.01_gr bl 50k 600s large 0d rs 2200_' ##part of the long filename to be removed. needed if peaklists are wanted in seperate files. Instead of the long name they will just be acq_x_peaks.txt

"""for baseline correction """
l=10^2
p=0.001 ##I use l=10^2 and p=0.001 here which seem to work fine for your data. have a play in the Raman_plot_singlefile programm to tune λ and p and see what they do


"""for peak picking"""
minimum_int=2000 ##specify here the minimum intensity for a peak to be considered a peak


""" for data plotting"""
colormap='rainbow' #colormap for 3D plot. order of colors can be reversed by adding _r after the name of the map. check here for available maps: https://matplotlib.org/examples/color/colormaps_reference.html

##formatting the graph
fig = plt.figure(figsize=(15,10))
ax = fig.gca(projection='3d')
ax.view_init(20, -80) ## change the orientation at which the 3D plot is viewed

""" if annotating of peaks isn't desired mark line 118-120 and press ctrl+1 to "deactivate" the lines by making them comments """
font_title=20 ##fontsize for axes titels
font_data=12 ##fontsize for data labels in graph
experiment='test' ##name under which the graph will be saved later
picture_format='.png' ##supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff





##defining a program to substract the baseline from our data. 
##I just found this online so don't ask me what exactly it does ;P I can only tell you it fits a polynomial to your data with asymmetric least squares
## source: https://stackoverflow.com/questions/29156532/python-baseline-correction-library
def baseline_als_optimized(y, lam, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z



##showing peaks annotated in 3D plot. again something I found online searching stackoverflow 
##source: https://stackoverflow.com/a/42915422/11777633    
class Annotation3D(Annotation):
    '''Annotate the point xyz with text s'''

    def __init__(self, s, xyz, *args, **kwargs):
        Annotation.__init__(self,s, xy=(0,0), *args, **kwargs)
        self._verts3d = xyz        

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.xy=(xs,ys)
        Annotation.draw(self, renderer)
        
def annotate3D(ax, s, *args, **kwargs):
    '''add anotation text s to to Axes3d ax'''

    tag = Annotation3D(s, *args, **kwargs)
    ax.add_artist(tag)


##reading in your data and substracting the baseline also finding peaks
for i,filename in enumerate(data_files):
    wave,intensity=np.loadtxt(filename, delimiter='\t',comments='#').T[[0,1]]
    baseline=baseline_als_optimized(intensity,l,p,niter=10) 
    int_corrected=intensity-baseline
    
    """peak picking"""
    name=os.path.basename(filename)
    name,ext=os.path.splitext(name)
    peaks,intense=signal.find_peaks(int_corrected,minimum_int)
    peaklist=np.array([wave[peaks],int_corrected[peaks]]).T
    np.savetxt(data_dir+re.sub(series, '',name)+'_peaks.txt',peaklist,delimiter='\t')
    xyzp=zip(wave[peaks],wave[peaks]*0+i/1,int_corrected[peaks])##needed for annotating data in plot
    

    """ plotting your data"""
    color =  matplotlib.cm.get_cmap(colormap)(i/data_file_count) 
    plt.plot(wave,wave*0+i/1,int_corrected,color=color)##0+i/1 = first spectrum at the front, 0-i/1 = first spectrum at the back

    for i,xyz_ in enumerate (xyzp):
      annotate3D(ax, s=int(wave[peaks][i]), xyz=xyz_, fontsize=font_data, xytext=(-3,3),
               textcoords='offset points', ha='left',va='bottom')      
    
    ##plot formatting
    ax.set_yticklabels([])
    ax.tick_params(axis='z', pad=20)
    ax.set_xlabel('Wavenumber $[cm^{-1}]$', fontsize=font_title, labelpad=20)
    ax.set_zlabel('Intensity [a.u.]', fontsize=font_title, labelpad=43)
    plt.savefig(data_dir+experiment+picture_format, transparent= True )

    

