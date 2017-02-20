# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 17:31:54 2017

@author: benk2002
"""
from __future__ import division

from pylab import *
from scipy import *
from scipy.signal import *
from numpy import *
from scipy import constants as k
import matplotlib.cm as cmaps
import sys
sys.path.append('c:/codes/pyHegel')
from pyHegel.util import readfile


class CplotV(object):
     """n numero du sweep
        mv nombre de points en V 
        s = 1 for save"""
    def __init__(self, filename, n,mv, fi, fm, s):
        self.filename = filename
        self.n = n        
        self.fi = fi
        self.fm = fm
        self.nc = 2**17
        self.s = s     
        self.x = np.arange(n*self.mv,(n+1)*self.mv)
        self.xs = shape(self.x)
        self.CMat = np.zeros((self.xs[0],self.nc,2), dtype=complex)
        
        
class CplotB(object):
    def __init__(self, filename, n, mb, mv, fi, fm, s):
        """n numero du sweep
        mb nombre de points en B
        mv nombre de points en V 
        s = 1 for save"""
        self.filename = filename
        self.n = n        
        self.fi = fi
        self.fm = fm
        self.nc = 2**17
        self.s = s    
        self.mb = mb
        self.mv = mv
        self.x =  np.arange(n,self.mb*self.mb,self.mv)
        self.xs = shape(self.x)
        self.CMat = np.zeros((self.xs[0],self.nc,2), dtype=complex)  
        
        
        
    def loadM(self):
        fc = lambda s: complex(s.replace('+-', '-'))
        for i in range (self.xs[0]):
            self.CMat[i]=loadtxt(self.filename.format(self.x[i]),converters={0:fc, 1:fc}, dtype=complex)
#            CBHP1-NoiseVDCB-V119-40db-I9-1M-BaseT_20170206-182343_d3_{:04.0f}.txt


        f = real(self.CMat[0,:,0])
        I = abs(sqrt(self.CMat[:,0,1])*1e-6/self.nc)
        S = self.CMat[:,:,1]

        self.fp = abs(f[0:self.fm])
        self.Sp = log(S[:,0:self.fm])
        self.Ip = I*1e9
        self.Ip[0:199] = -self.Ip[0:199]
        
        if self.s == 1:
            Mat=(f,I,S)
            save('CMatN{:02.0f}'.format(self.n),Mat)
            
       
        
    def plotM(self):
        fig = figure(figsize = [16,9])
        ax1 = fig.add_subplot(1,1,1)
        plt.set_cmap(cmaps.viridis)

        for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
            ax1.get_xticklabels() + ax1.get_yticklabels()):
            item.set_fontsize(20)

#        CS1 = ax1.contourf(self.fp , self.Ip , self.Sp , self.l , vmin = -8.046, vmax = -2)
        CS1 = ax1.pcolormesh(self.fp , self.Ip , self.Sp , vmin = -8.046, vmax = -2)


            #cbar=fig.colorbar(CS1, ax=ax1, shrink=0.9)
            #ax1.set_xlim(xmin = 4.173, xmax = 4.88)
        plt.xlabel ("$f(Hz)$")
        plt.ylabel (r"$I_{bias}(nA)$")
        #cbar.ax.set_ylabel ("$\sigma (e^2/h)$ ")
        plt.tight_layout()

        fig.savefig("CPNoiseN{:02.0f}.jpg".format(self.n))
        plt.clf()
        plt.close()
        del CS1 , ax1 , self.fp , self.Ip , self.Sp , self.l

