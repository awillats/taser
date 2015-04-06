# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:47:31 2015

@author: Adam
"""


"""
set up defaults!
probe block setParams
help text

updateparams function improve

better default sigBlock
"""
import numpy as np
from numpy import random as npr
import matplotlib.pyplot as plt
from pylab import *
import scipy.signal as signal

class SigBlock:
    SAMPLE_FREQ = float(25E3)
    
    def __init__(self, sigType, tStart, tLen, params):
        self.sigType = sigType
        self.tStart = tStart
        self.tLen = tLen
        self.params = params
        
    def genArraySetup(self):
        self.nLen = self.tToN(self.tLen)

    def genArray(self):
        self.genArraySetup()
        self.valArray = np.zeros(self.nLen)
    
    def getTArray(self):
        return np.linspace(self.tStart,self.tStart+self.tLen,num=self.tToN(self.tLen)) 
   
    def plot(self):
        plt.plot(self.getTArray(),self.valArray)

    def scaleShiftArray(self,scaleVal=1,shiftVal=0):
        if not(hasattr(self,'valArray')):
            print "value array wasn't generated, but is zeros"
            self.genArray()
            
        self.valArray *= scaleVal
        self.valArray += shiftVal
        
    def scaleShiftArraySection(self,secStartT,sectionLenT,scaleVal=1,shiftVal=0):
        if not(hasattr(self,'valArray')):
            print "value array wasn't generated, but is zeros"
            self.genArray()
        startN = self.tToI(secStartT)       
        endN = self.tToI(secStartT+sectionLenT)
       # print 'shifting from',secStartT,' to',secStartT+sectionLenT 
       # print startN,endN, endN-startN, len(shiftVal), len(self.valArray)
        
        self.valArray[startN:endN] *= scaleVal
        self.valArray[startN:endN] += shiftVal   
        
        
    def boundArray(self,lowBound,highBound):
        np.clip(self.valArray,lowBound,highBound,out=self.valArray)
        return self.valArray
        
    def HPLPArray(self,filtParams):
        #http://mpastell.com/2010/01/18/fir-with-scipy/
        filtN = filtParams['n'] if filtParams.has_key('n') else 100
        win = filtParams['winType'] if filtParams.has_key('winTpye') else 'blackmanharris'
        fc=0
        if filtParams.has_key('LPcutHz'):
            if filtParams.has_key('HPcutHz'):
                fc = BPcoeffs(filtN, filtParams['LPcutHz']/self.SAMPLE_FREQ,
                              filtParams['HPcutHz']/self.SAMPLE_FREQ, win)
            else:
                fc = LPcoeffs(filtN, filtParams['LPcutHz']/self.SAMPLE_FREQ, win)
        else:
            if filtParams.has_key('HPcutHz'):
                fc = HPcoeffs(filtN, filtParams['HPcutHz']/self.SAMPLE_FREQ, win)
            else:
                print " 'what's going on' - marvin gaye. No filter params"
                return self.valArray
        
        self.valArray = filtSig(self.valArray,fc)
        return self.valArray    
        
    def endT(self):
        return self.tStart+self.tLen

    
    def tToN(self,t):
        nF = float(self.SAMPLE_FREQ*t);
        n = int(np.round(nF))
        
        if (n-nF > 1e-10):
            print "note, the segment length is not an integer number of samples"
            print "using t =", n/self.SAMPLE_FREQ, "not", nF/self.SAMPLE_FREQ
            print 'ns:',n,nF
            print 'diff:',n-nF
            
        return n
        
    def tToI(self,t):
        return self.tToN(t-self.tStart)

    def nyquistCheck(self,f):
        if f > 2*self.SAMPLE_FREQ:
            print 'warning, ', f,' is greater than the nyquist rate:', 
            2*self.SAMPLE_FREQ, 'Hz'
            
    def printToTxt(self,fName,doSaveMeta=False):
        fName1 = fName+'.txt'
        #fmt string
        np.savetxt(fName1,self.valArray)
        
        if doSaveMeta:
            fNameMeta = fName+'_meta.txt'
            f = open(fNameMeta,'w')
            m = str(self.params)
            f.write(m)
            f.close()
        
        
        
def chooseBlock(blockType, params):
    print blockType
    blocks = {'WN': WNBlock,

                'Ramp': RampBlock,

                'Pulse': PulseBlock,

                'Sine': SineBlock,

                'Probe': ProbeBlock
     }
     
    return blocks[blockType](params)     
     
class WNBlock(SigBlock):
    def __init__(self,p):
        self.setFromParams(p)
        SigBlock.__init__(self,'WN',self.tStart,self.tLen,p)
        
    def setFromParams(self,p):
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.mean = p['mean']
        self.var = p['var']
    def genArray(self):

        if self.sigType == 'Static':
            print 'skipping gen, because this block is static'
            return
        
        SigBlock.genArraySetup(self)

        if hasattr(self,'wnSeed'):
            #print 'seed is', self.wnSeed
            npr.seed(self.wnSeed)
        else:
            #print 'no seed'
            npr.seed(None)
            
        self.valArray = npr.normal(loc=self.mean,scale=sqrt(self.var),size=self.nLen)
        return self.valArray
             
        
    def setRandSeed(self,seed):
        self.wnSeed = seed
        self.params.update({'seed':seed})

class RampBlock(SigBlock):
    def __init__(self,p):
        self.setFromParams(p)
        SigBlock.__init__(self,'Ramp',self.tStart,self.tLen,p)
    def setFromParams(self,p):
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.amp = p['amp']
        self.dir = p['dir']
    def genArray(self):
        SigBlock.genArraySetup(self)

        if self.dir.upper()=='UPRIGHT':
            self.valArray = np.linspace(0,self.amp,num=self.nLen)
        else:
            self.valArray = np.linspace(self.amp,0,num=self.nLen)
   
        return self.valArray

    def setParams(self, rampH, rampDir):
        self.amp = rampH
        self.dir = rampDir
        self.params.update({'dir':rampDir,'amp':rampH})
        
class PulseBlock(SigBlock):
    def __init__(self,p):
        self.setFromParams(p)
        SigBlock.__init__(self,'Pulse',self.tStart,self.tLen,p)
    def setFromParams(self,p):

        if (not p.has_key('lowT')) and p.has_key('tLen'):
            p['hiT'] = p['tLen']
            p['lowT'] = 0
            
        
        tLen = p['hiT']+p['lowT']
        self.tLen = tLen
        p['tLen'] = tLen
        self.tStart = p['tStart']
        self.hiT = p['hiT']
        self.lowT = p['lowT'] #=0
        self.hiVal = p['hiVal']
    def genArray(self):
        SigBlock.genArraySetup(self)
            
        hiN = self.tToN(self.hiT)
        lowN = self.tToN(self.lowT)
        
        self.valArray = np.concatenate((self.hiVal*np.ones(hiN),np.zeros(lowN)))
        return self.valArray
              
    def setParams(self, hiT,lowT, hiVal):
        self.hiT = hiT
        self.lowT = lowT
        self.hiVal = hiVal
        self.params.update({'hiT':hiT,'lowT':lowT, 'hiVal':hiVal})                  

class SineBlock(SigBlock):
    def __init__(self,p):
        self.setFromParams(p)
        SigBlock.__init__(self,'Sine',self.tStart,self.tLen,p)
        self.nyquistCheck(p['freqHz'])

    def setFromParams(self,p):
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.amp = p['amp']
        self.freq = p['freqHz']
        self.dcOff = p['dcOff']
        self.phase0 = p['phase0']
    def genArray(self):
        SigBlock.genArraySetup(self)
        phaseEnd = self.phase0+self.tLen*self.freq*2*np.pi
        phaseArray = np.linspace(self.phase0,phaseEnd,num=self.nLen)
        self.valArray = self.amp*np.sin(phaseArray)+self.dcOff
        return self.valArray
              
    def setParams(self, amp,freqHz,dcOff=0,phase0=0):
        self.amp = amp
        self.freq = freqHz
        self.dcOff = dcOff
        self.phase0 = phase0
        self.params.update( {'amp':amp,'freq':freqHz,'dcOff':dcOff,'phase0':phase0})                  

class ProbeBlock(SigBlock):
    """ Gaussian probe stimulus """
    def __init__(self,p):
        self.setFromParams(p)
        SigBlock.__init__(self,'Probe',self.tStart,self.tLen,p)
    def setFromParams(self,p):
        self.tLen = p['tWide']
        self.tStart = p['tCenter'] - self.tLen/2.0
        p['tStart'] = self.tStart
        p['tLen'] = self.tLen
        print "block starts at tCenter-tLen/2 > tStart=",self.tStart
        self.amp = p['amp']
        self.std = p['std']
        self.tCenter = p['tCenter']
    def genArray(self):
        SigBlock.genArraySetup(self)
        t = self.getTArray()
        
        self.valArray = self.amp*gaussian(t,self.tCenter,self.std)
        return self.valArray   
        
    def setParams(self):
        print 'fix me'
        self.params.update( {'amp':self.amp})

#just math    
def gaussian(x, mu, sig):
   return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def LPcoeffs(n, cutoffRatio,win='hanning'):
    a = signal.firwin(n,cutoff=cutoffRatio,window=win)
    return a
def HPcoeffs(n, cutoffRatio,win='hanning'):
    a = signal.firwin(n,cutoff=cutoffRatio,window=win)
    a = -a
    a[n/2]= a[n/2]+1
    return a
def BPcoeffs(n,lpCutRatio,hpCutRatio,win='blackmanharris'):
    a = LPcoeffs(n,lpCutRatio,win)
    b = HPcoeffs(n,hpCutRatio,win)
    d = - (a+b); d[n/2] = d[n/2] + 1
    return d
def filtSig(sig,filt):
    return signal.lfilter(filt,1.0,sig)   
    
def mfreqz(b,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * log10 (abs(h))
    subplot(211)
    plot(w/max(w),h_dB)
    ylim(-150, 5)
    ylabel('Magnitude (db)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Frequency response')
    subplot(212)
    h_Phase = unwrap(arctan2(imag(h),real(h)))
    plot(w/max(w),h_Phase)
    ylabel('Phase (radians)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Phase response')
    subplots_adjust(hspace=0.5)

#Plot step and impulse response
def impz(b,a=1):
    l = len(b)
    impulse = repeat(0.,l); impulse[0] =1.
    x = arange(0,l)
    response = signal.lfilter(b,a,impulse)
    subplot(211)
    stem(x, response)
    ylabel('Amplitude')
    xlabel(r'n (samples)')
    title(r'Impulse response')
    subplot(212)
    step = cumsum(response)
    stem(x, step)
    ylabel('Amplitude')
    xlabel(r'n (samples)')
    title(r'Step response')
    subplots_adjust(hspace=0.5)