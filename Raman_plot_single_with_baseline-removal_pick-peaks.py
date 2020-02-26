# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:57:15 2020

@author: ms01106
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy import signal 

data_dir = "//surrey.ac.uk/personal/HS129/ms01106/.System/Desktop/Katrin Raman Scripts/" 
wave,intensity=np.loadtxt(data_dir+'2020.01_gr bl 50k 600s large 0d rs 2200_Acq28.txt', delimiter='\t',comments='#').T[[0,1]]


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

"""There are two parameters: p for asymmetry and λ for smoothness. 
Both have to be tuned to the data at hand. 
We found that generally 0.001 ≤ p ≤ 0.1 is a good choice (for a signal with positive peaks) and 10^2 ≤ λ ≤ 10^9 , but exceptions may occur. 
In any case one should vary λ on a grid that is approximately linear for log λ
NOTE: There seems to be a paper about this by P.Eilers and H.Boelens from 2005 that explains the algorithm more but I couldn't find it.... """


baseline=baseline_als_optimized(intensity,10^2,0.0001,niter=10)
int_corrected=intensity-baseline

peaks,intense=signal.find_peaks(int_corrected,5000)
print(wave[peaks])


plt.plot(wave,intensity,'r--',label='original data')
plt.plot(wave,baseline,'g',label='fitted baseline')
plt.legend(fontsize=12)
plt.show()

plt.plot(wave,int_corrected,'b',label='baseline-subtracted data')
plt.plot(wave[peaks],int_corrected[peaks],'kx',label='selected peaks')
plt.legend(fontsize=12)
plt.plot()
