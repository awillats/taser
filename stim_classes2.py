# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:47:31 2015

@author: Adam
"""


"""
eliminate extra deepcopies


probe from center, from start

~help text
~updateparams function improve
Xset up defaults!
Xbetter default sigBlock
Xprobe block setParams
Xplot trigs->in stream

updateFromP needs deepcopy?
pulse  and probe NEED certain params, don't default well
?self.tStart => p = self.params. p['tStart']

"""

import numpy as np
from numpy import random as npr
import matplotlib.pyplot as plt

#for filtering
from pylab import *
import scipy.signal as signal

import copy
import shutil
import os

class SigBlock:
    """
    The building block for a stimulus, by default all 0s
    
    ...
    
    Attributes
    ----------
    SAMPLE_FREQ : float
        sample frequency [Hz]
    sigType : String
        - specifies the type of block
        - '' for silence
        - 'Static' for a compiled static signal, not dynamically regenerated
        - 'WN' for white noise
        - 'Pulse' for an on/off pulse
        - 'Sine' for a sine wave
        - 'Probe' for a gaussian waveform
    params : dict
        tStart : float
            start time [seconds]
        tLen : float
            length of block [seconds]
        val : float
            by default assumes the signal=val for the whole time
        
    
    Methods
    ----------
    modTheseParams()
        newP
    getTArray()
    
    genArray()
    
    plot()
    
    scaleShiftArray()
        scaleVal=1
        shiftVal=0
        
    scaleShiftArraySection()
        secStartT
        sectionLenT
        scaleVal=1
        shiftVal=0
        
    boundArray()
        lowBound,highBound
        
    HPLPArray()
        {'n','winType','LPcutHz','HPcutHz'}
    endT()
    
    tToN()
        t
    tToI()
        t
    nyquistCheck()
        f
    printToTxt()
        fNamePre
        doSaveMeta=false
    
    
    """
    SAMPLE_FREQ = float(25E3)
        
    DEF_P = {'tStart':0,'tLen':1,'val':0}
    
    def __init__(self, sigType='',params={}):
        self.sigType = sigType
        #set params to default
        self.params = copy.deepcopy(self.DEF_P)
        #override defaults with incoming values
        self.modTheseParams(params)
        
        
    def genArraySetup(self):
        self.nLen = self.tToN(self.tLen)
        
   
    def modTheseParams(self,newP_):
        ''' updates only the parameters entered '''
        newP = copy.deepcopy(newP_)
        for k,v in newP.items():
            self.params[k] = v 
        self._updateFromP(self.params)
    def shiftParams(self,shiftP_):
        #FIX ME NOT NECESSARy     
        self._updateFromP(copy.deepcopy(self.params))        
        ''' updates only the parameters entered '''
        shiftP = copy.deepcopy(shiftP_)
        for k,v in shiftP.items():
            self.params[k] += v 
            
        self._updateFromP(self.params)  
    def shiftStart(self,tShift):
        ssp = {'tStart':tShift}
        self.shiftParams(copy.deepcopy(ssp))
    def replaceParams(self,newP):
        '''deletes existing params, replaces only with input params '''
        print 'why'
        self.params = copy.deepcopy(newP)
        self._updateFromP(self.params)   
        
    def _updateFromP(self,p):
        ''' sets shortcut variables for params '''
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.val = p['val']
        
    def genArray(self):
        ''' generate values '''
        self.genArraySetup()
        self.valArray = np.zeros(self.nLen)+self.val
        #+self.params['val']
    def _testGenValArray(self):
        if not(hasattr(self,'valArray')):
            print "value array wasn't generated, but is now"
            self.genArray()
            
    def getTArray(self):
        '''prints out an array from tStart to tStart+tLen using the sample freq'''
        return np.linspace(self.tStart,self.tStart+self.tLen,num=self.tToN(self.tLen)) 
   
    def plot(self):
        '''plots valArray'''
        self._testGenValArray()
        print 'lens',len(self.getTArray()), len(self.valArray)
        plt.plot(self.getTArray(),self.valArray)
   
    def scaleShiftArray(self,scaleVal=1,shiftVal=0):
        '''takes whole valArray, multiplies by scaleVal, THEN adds shiftVal'''
        self._testGenArray()
            
        self.valArray *= scaleVal
        self.valArray += shiftVal
        
    def scaleShiftArraySection(self,sectStartT,sectLenT,scaleVal=1,shiftVal=0):
        self._testGenValArray()
        
        startN = self.tToI(sectStartT)       
        endN = self.tToI(sectStartT+sectLenT)
        #print 'shifting from',sectStartT,' to',sectStartT+sectLenT 
        #print startN,endN, endN-startN, len(shiftVal), len(self.valArray)
        
        self.valArray[startN:endN] *= scaleVal
        self.valArray[startN:endN] += shiftVal   
        
    def boundArray(self,lowBound,highBound):
        self._testGenValArray()
        np.clip(self.valArray,lowBound,highBound,out=self.valArray)
        return self.valArray
        
    def HPLPArray(self,filtParams):
        """
        asd
        Arguments
        --------
        n
        winType
        LPcutHz
        HPcutHz
        """
        #http://mpastell.com/2010/01/18/fir-with-scipy/
        filtN = filtParams['n'] if filtParams.has_key('n') else 1250
        win = filtParams['winType'] if filtParams.has_key('winTpye') else 'blackmanharris'
        
        fc = 0
        
        if filtParams.has_key('LPcutHz'):
            if filtParams.has_key('HPcutHz'):
                fc = BPcoeffs(filtN, filtParams['LPcutHz']/self.SAMPLE_FREQ,
                              filtParams['HPcutHz']/self.SAMPLE_FREQ, win)
            else:
                fc = LPcoeffs(filtN, filtParams['LPcutHz']/self.SAMPLE_FREQ, win)
                print 'LP filter'
        else:
            if filtParams.has_key('HPcutHz'):
                fc = HPcoeffs(filtN, filtParams['HPcutHz']/self.SAMPLE_FREQ, win)
            else:
                #print " 'what's going on' - marvin gaye. No filter params"
                return self.valArray
        
        self.valArray = filtSig(self.valArray,fc)
        return self.valArray    
        
    def endT(self):
        '''prints end point'''
        return self.tStart+self.tLen

    
    def tToN(self,t):
        ''' converts a length of time to a number of samples using the sample freq.'''
        nF = float(self.SAMPLE_FREQ*t);
        n = int(np.round(nF))
        
        if (n-nF > 1e-10):
            print "note, the segment length is not an integer number of samples"
            print "using t =", n/self.SAMPLE_FREQ, "not", nF/self.SAMPLE_FREQ
            print 'ns:',n,nF
            print 'diff:',n-nF
            
        return n
        
    def tToI(self,t):
        '''converts a time to a number of samples since the start of the block using the sample freq.'''
        return self.tToN(t-self.tStart)

    def nyquistCheck(self,f):
        ''' checks if a frequency is > f_sample/2 '''
        if f > self.SAMPLE_FREQ/2.0:
            print 'warning, ', f,' is greater than the nyquist freq.:', 
            self.SAMPLE_FREQ/2.0, 'Hz'
    def assignFolder(self,foldName):
        self.foldName = foldName
        makeFolder(foldName)
        
        
    def printToTxt(self,fNamePre,doSaveMeta=False):
        if hasattr(self,'foldName'):
            fNamePre = self.foldName+fNamePre
        
        fName = fNamePre+'.txt'
        #fmt string
        np.savetxt(fName,self.valArray)
        
        if doSaveMeta:
            fNameMeta = fNamePre+'_meta.txt'
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
    """
    Block of gaussian white noise
    ...
    
    Params
    ----------

    p : dict 
        tStart : float
            block start [seconds]
            
        tLen : float
             length of block [seconds]
             
        mean : float
             mean
             
        var : float
             variance
    Methods
    --------
    setRandSeed():
        
    Defaults
    ----------
    {'tStart':0,'tLen':1,'mean':0,'var':1,'seed':None}
    """
    DEF_P = {'tStart':0,'tLen':1,'mean':0,'var':1,'seed':None}
    
    def __init__(self,p={}):
        SigBlock.__init__(self,'WN',p)
        
    def _updateFromP(self,p):
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.mean = p['mean']
        self.var = p['var']
        self.seed = p['seed']
        
    def genArray(self):
        self.genArraySetup()
        
        #if self.seed == None: print 'no seed'
        npr.seed(self.seed)
        ###self.valArray = npr.normal(loc=self.mean,scale=sqrt(self.var),size=self.nLen)
        self.valArray = npr.normal(loc=self.mean,scale=sqrt(self.var),size=self.nLen)
        
        if self.params.has_key('lBound'):
            self.boundArray(self.params['lBound'],self.params['uBound'])
        return self.valArray
 
    def setRandSeed(self,seed):
        ''' sets seed for random number generation '''
        self.modTheseParams({'seed':seed})

class RampBlock(SigBlock):
    """
    sloped rise/fall
    ...
    
    Params
    ----------

    p : dict 
        tStart : float
            block start [seconds]
            
        tLen : float
             length of block [seconds]
             
        amp : float
             max value of the ramp
             
        dir : String
            - if dir == 'UPRIGHT' the ramp slopes upwards from 0 to amp
            - if dir is not = 'UPRIGHT' the ramp slopes downwards from amp to 0
    """
    DEF_P = {'tStart':0,'tLen':1,'amp':1,'dir':'UPRIGHT'}
    def __init__(self,p={}):
        SigBlock.__init__(self,'WN',p)
        
    def _updateFromP(self,p):
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.amp = p['amp']
        self.dir = p['dir']
    def genArray(self):
        self.genArraySetup()

        if self.dir.upper()=='UPRIGHT':
            self.valArray = np.linspace(0,self.amp,num=self.nLen)
        else:
            self.valArray = np.linspace(self.amp,0,num=self.nLen)
   
        return self.valArray
     
class PulseBlock(SigBlock):
    DEF_P = {'tStart':0,'hiT':1,'lowT':0,'hiVal':1}
    
    def __init__(self,p={}):
        """
        if not p.has_key('hiT') and not p.has_key('lowT') and not p.has_key('tLen'):
            p['hiT'] = self.DEF_P['hiT']
            p['lowT'] = self.DEF_P['lowT']
            """
        #FIX ME
        
        SigBlock.__init__(self,'Pulse',p)
        #set as defaults, mod parameters given as inputs, 
    def _hiLowLen(self,p):
        try:
            #why doesn't try work
            if p.has_key('tLen'):
                if p.has_key('hiT') and not p.has_key('lowT'):
                    p['lowT'] = p['tLen']-p['hiT']
                elif p.has_key('lowT') and not p.has_key('hiT'):
                    p['hiT'] = p['tLen']-p['lowT']
                elif not p.has_key('hiT') and not p.has_key('lowT'):
                    p['hiT'] = p['tLen']
                    p['lowT'] = 0
            else:
                if not p.has_key('hiT') and not p.has_key('lowT'):
                    p['hiT'] = self.DEF_P['hiT']
                    p['lowT'] = self.DEF_P['lowT']
                    p['tLen'] = self.DEF_P['lowT'] + self.DEF_P['hiT']
                if not p.has_key('hiT'):
                    p['hiT'] = self.DEF_P['hiT']
                if not p.has_key('lowT'):
                    p['lowT'] = self.DEF_P['lowT']
                p['tLen'] = p['hiT']+p['lowT']
            return p
        except e:
            #print 'not enough params', e
            return {}
            
    def modTheseParams(self,p):
        '''overwrites from SigBlock'''
        #p = self._hiLowLen(p)
        SigBlock.modTheseParams(self,p) 
    def shiftParams(self,p):
        '''overwrites from SigBlock'''
        #p = self._hiLowLen(p)
        SigBlock.shiftParams(self,p)
                
    def _updateFromP(self,p):
        p=self._hiLowLen(p)
        self.tLen = p['tLen']        
        self.tStart = p['tStart']
        self.hiT = p['hiT']
        self.lowT = p['lowT']
        self.hiVal = p['hiVal']
    def genArray(self):
        self.genArraySetup()
            
        hiN = self.tToN(self.hiT)
        lowN = self.tToN(self.lowT)
        
        self.valArray = np.concatenate((self.hiVal*np.ones(hiN),np.zeros(lowN)))
        return self.valArray
                             

class SineBlock(SigBlock):
    
    DEF_P = {'tStart':0,'tLen':1,'amp':1,'freqHz':10,'dcOff':0,'phase0':0}    
    
    def __init__(self,p={}):
        SigBlock.__init__(self,'Sine',p)
        self.nyquistCheck(p['freqHz'])

    def _updateFromP(self,p):
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.amp = p['amp']
        self.freq = p['freqHz']
        self.dcOff = p['dcOff']
        self.phase0 = p['phase0']
    def genArray(self):
        self.genArraySetup()
        phaseEnd = self.phase0+self.tLen*self.freq*2*np.pi
        phaseArray = np.linspace(self.phase0,phaseEnd,num=self.nLen)
        self.valArray = self.amp*np.sin(phaseArray)+self.dcOff
        return self.valArray
                         

class ProbeBlock(SigBlock):
    """ Gaussian probe stimulus """
    DEF_P = {'tCenter':0,'tWide':1,'amp':1,'std':1.0/6.0}  
    def __init__(self,p={}):
        #self._centerToStart(p)
        SigBlock.__init__(self,'Probe',p)
    def _centerToStart(self,p):
        p['tLen'] = p['tWide']
        p['tStart'] = p['tCenter'] - p['tLen']/2.0
        print "block starts at tCenter-tLen/2 > tStart=",p['tStart']
        return p
    def shiftStart(self,tShift):
        self.shiftParams({'tCenter':tShift})   
    def _updateFromP(self,p):
        self._centerToStart(p)
        self.tStart = p['tStart']
        self.tLen = p['tLen']
        self.amp = p['amp']
        self.std = p['std']
        self.tCenter = p['tCenter']
    def genArray(self):
        #gBlock.genArraySetup(self)
        self.genArraySetup()
        t = self.getTArray()
        self.valArray = self.amp*gaussian(t,self.tCenter,self.std)
        return self.valArray   
        

#just math
def makeFolder(foldName,doReplace=False):
    if not os.path.exists(foldName):
        os.makedirs(foldName)
    else: 
        if doReplace:
            print 'dir rewritten'
            shutil.rmtree(foldName)
            os.makedirs(foldName)
        
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
    a =  LPcoeffs(n,lpCutRatio,win)
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