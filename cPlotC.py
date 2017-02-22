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
import re
from scipy import constants as k
import matplotlib.cm as cmaps
import time
import sys
sys.path.append('c:/codes/pyHegel')
from pyHegel.util import readfile



        
        
class Cplot(object):
    def __init__(self, filename, n, vb, fi, fm, s):
        """n numero du sweep
        mb nombre de points en B
        mv nombre de points en V 
        s = 1 for save"""
        self.filename = filename
        self.fdname = filename[:-4] + '_d3_{:04.0f}.txt'
        self.n = n        
        self.fi = fi
        self.fm = fm
        self.nc = 2**17
        self.s = s  
        self.vb = vb
        self.V=0
        self.B=0
        
        
    
    
        
    def loadR(self):
        fc = lambda s: complex(s.replace('+-', '-'))
        VB = readfile(self.filename, getheaders = True)
        st = VB[2][5]
        a = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*",st)
        mb = int(a[0])
        mv = int(a[1])
        M = VB[0]
        
        if shape(M)[1] != mb:
            l = int(shape(M[1])[0]/ mv)
            x = np.arange(l*mv,shape(M[1])[0])            
            M = np.delete(M,x, axis=1)            
            M = np.reshape(M,(shape(M)[0], l, mv))
            
        
            
        if self.vb == 0:
            x =  np.arange(self.n, mb*mv, mv)
            self.B = M[0]
            self.V = M[1]
            
        else:
            x = np.arange(self.n*mv,(self.n+1)*mv)
            self.B = M[1]
            self.V = M[0]
             
        xs = shape(x)
            
        self.CMat = np.zeros((xs[0],self.nc,2), dtype=complex)  
        
        for i in range (xs[0]):
            self.CMat[i]=loadtxt(self.fdname.format(x[i]),converters={0:fc, 1:fc}, dtype=complex)
#            print('loading:' i*100/xs[0] '%')
            
            update_progress(i/xs[0])
#            CBHP1-NoiseVDCB-V119-40db-I9-1M-BaseT_20170206-182343_d3_{:04.0f}.txt


        f = real(self.CMat[0,:,0])
        I = abs(sqrt(self.CMat[:,0,1])*1e-6/self.nc)
        S = self.CMat[:,:,1]

        self.fp = abs(f[0:self.fm])
        self.Sp = log(S[:,0:self.fm])
        self.Ip = I*1e9
        self.Ip[0:199] = -self.Ip[0:199]
        
        self.Mat=(S)
        self.MatP=(f, self.B, self.V, I)
        
        if self.s == 1:  
            if self.vb == 1:
                save(self.filename[:-4] + 'Vs_B={:02.3f}T'.format(self.B[0,self.n]),self.Mat)
                save(self.filename[:-4] + 'Vs_B={:02.3f}T'.format(self.B[0,self.n]),self.MatP)
            else:      
                save(self.filename[:-4] + 'Bs_V={:02.3f}V'.format(self.V[0,self.n]),self.Mat)
                save(self.filename[:-4] + 'Bs_V={:02.3f}V'.format(self.V[0,self.n]),self.MatP)
                
            
    def loadM(self):
        self.Mat=load(self.filename)
        
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
        
        
    def Stat(self):
        """ compute the noise power spectrum"""
        
        for i in range (shape(self.Mat[4])[0]): 
           
            X = np.delete(self.Mat[4][i],np.where(abs(f)<11000))
            fX = np.delete(self.Mat[0],np.where(abs(f)<11000))
    
            X = np.delete(X,np.where(abs(fX)>62000))
            fX = np.delete(fX,np.where(abs(fX)>62000))
    
            X = np.delete(X,np.where((abs(fX)>14300) & (abs(fX)<20600)))
            fX = np.delete(fX,np.where((abs(fX)>14300) & (abs(fX)<20600)))
    
            X = np.delete(X,np.where((abs(fX)>20810) & (abs(fX)<20890)))
            fX = np.delete(fX,np.where((abs(fX)>20810) & (abs(fX)<20890)))
    
            X = np.delete(X,np.where((abs(fX)>36700) & (abs(fX)<37300)))
            fX = np.delete(fX,np.where((abs(fX)>36700) & (abs(fX)<37300)))
    
            X = np.delete(X,np.where((abs(fX)>39200) & (abs(fX)<40600)))
            fX = np.delete(fX,np.where((abs(fX)>39200) & (abs(fX)<40600)))
    
            X = np.delete(X,np.where((abs(fX)>42200) & (abs(fX)<43500)))
            fX = np.delete(fX,np.where((abs(fX)>42200) & (abs(fX)<43500)))

            M2n[i] = np.sum(abs(X))
            SX[i] = shape(X)[0]
            
        MatN = (self.Mat[1], SX, M2n)
        return MatN
 



def update_progress(progress):
        barLength = 10 # Modify this to change the length of the progress bar
        status = ""
        if isinstance(progress, int):
            progress = float(progress)
        if not isinstance(progress, float):
            progress = 0
            status = "error: progress var must be float\r\n"
        if progress < 0:
            progress = 0
            status = "Halt...\r\n"
        if progress >= 1:
            progress = 1
            status = "Done...\r\n"
        block = int(round(barLength*progress))
        text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
        sys.stdout.write(text)
        sys.stdout.flush()
           
        
#class CplotV(CplotB):        
#    def __init__(self, filename, n, mv, fi, fm, s):
#        """n numero du sweep
#        mv nombre de points en V 
#        s = 1 for save"""
#        self.filename = filename
#        self.fdname = filename[:-4] + '_d3_{:04.0f}.txt'
#        self.n = n        
#        self.fi = fi
#        self.fm = fm
#        self.mv = mv
#        self.nc = 2**17
#        self.s = s     
#        self.x = np.arange(n*self.mv,(n+1)*self.mv)
#        self.xs = shape(self.x)
#        self.CMat = np.zeros((self.xs[0],self.nc,2), dtype=complex)

