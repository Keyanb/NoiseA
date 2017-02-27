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
import os
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
            M = swapaxes(M, 1, 2)           
            self.B = transpose(M[0])
            self.V = transpose(M[1])   
        else:
            self.B = M[0]
            self.V = M[1]
       
        if self.vb == 0:
            x =  np.arange(self.n, mb*mv, mv)
                       
        else:
            x = np.arange(self.n*mv,(self.n+1)*mv)
             
        xs = shape(x)
            
        self.CMat = np.zeros((xs[0],self.nc,2), dtype=complex)  
        
        if self.vb == 1:
            fna = self.filename[:-4] + 'Vs_B={:02.3f}T.npy'.format(self.B[self.n,0])
            fnaP = self.filename[:-4] + 'Vs_B={:02.3f}T-P.npy'.format(self.B[self.n,0])
        else:
            fna = self.filename[:-4] + 'Bs_V={:02.3f}V.npy'.format(self.V[0,self.n])
            fnaP = self.filename[:-4] + 'Bs_V={:02.3f}V-P.npy'.format(self.V[0,self.n])
        
        if os.path.isfile(fna):
            print('T')
            self.Mat = load(fna)
            self.MatP = load(fnaP)
            
        else:        
            for i in range (xs[0]):
                self.CMat[i]=loadtxt(self.fdname.format(x[i]),converters={0:fc, 1:fc}, dtype=complex)
                update_progress(i/xs[0])

            f = real(self.CMat[0,:,0])
            I = abs(sqrt(self.CMat[:,0,1])*1e-6/self.nc)
            S = self.CMat[:,:,1]
               
            self.Mat=(S)
            self.MatP=(f, self.B, self.V, I)
        
            if self.s == 1:  
                if self.vb == 1:
                    save(self.filename[:-4] + 'Vs_B={:02.3f}T'.format(self.B[self.n,0]),self.Mat)
                    save(self.filename[:-4] + 'Vs_B={:02.3f}T-P'.format(self.B[self.n,0]),self.MatP)
                else:      
                    save(self.filename[:-4] + 'Bs_V={:02.3f}V'.format(self.V[0,self.n]),self.Mat)
                    save(self.filename[:-4] + 'Bs_V={:02.3f}V-P'.format(self.V[0,self.n]),self.MatP)
                
            del f, I, S, x, xs, VB, M, fc, mb, mv, a, st
        
    def plotM(self):
        fig = figure(figsize = [16,9])
        ax1 = fig.add_subplot(1,1,1)
        plt.set_cmap(cmaps.viridis)
        
        fp = abs(self.MatP[0][0:self.fm])
        Sp = log(self.Mat[:,0:self.fm])
        Ip = self.MatP[3]*1e9
        Ip[0:199] = -Ip[0:199]

        for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
            ax1.get_xticklabels() + ax1.get_yticklabels()):
            item.set_fontsize(20)

#        CS1 = ax1.contourf(self.fp , self.Ip , self.Sp , self.l , vmin = -8.046, vmax = -2)
        if self.vb == 1:
            CS1 = ax1.pcolormesh(fp , self.V[0,:] , Sp , vmin = -8.046, vmax = -2)
        else:
            CS1 = ax1.pcolormesh(fp , self.B[:,0] , Sp , vmin = -8.046, vmax = -2)

            #cbar=fig.colorbar(CS1, ax=ax1, shrink=0.9)
            #ax1.set_xlim(xmin = 4.173, xmax = 4.88)
        plt.xlabel ("$f(Hz)$")
        
        #cbar.ax.set_ylabel ("$\sigma (e^2/h)$ ")
        plt.tight_layout()

        if self.vb == 0:
            plt.ylabel (r"$B(T)$")
            fig.savefig("CPNoiseBs_V={:02.3f}V.jpg".format(self.V[0,self.n]))
        else:
            plt.ylabel (r"$V_{bias}(V)$")
            fig.savefig("CPNoiseVs_B={:02.3f}T.jpg".format(self.B[self.n,0]))
        plt.clf()
        plt.close()
        del CS1 , ax1, fp, Sp, Ip 
        
        
    def Stat(self):
        """ compute the noise power spectrum"""
        
        f = self.MatP[0]
        fx = f
        M2n = np.zeros(shape(self.Mat)[0])
        SX = np.zeros(shape(self.Mat)[0])
        
        for i in range (shape(self.Mat)[0]): 
           
            X = np.delete(self.Mat[i],np.where(abs(f)<11000))
            fX = np.delete(fx,np.where(abs(f)<11000))
    
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
        
        plot(fX,X)
        MatN = (self.Mat[1], SX, M2n)
        if self.vb == 0:
            save("StatBs_V={:02.3f}V".format(self.V[0,self.n]), MatN)
        else:
            save("StatVs_B={:02.3f}T".format(self.B[0,self.n]), MatN)
        return MatN
        
        del X, fX
 



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

