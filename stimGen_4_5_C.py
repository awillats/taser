# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 15:32:03 2015

@author: Adam
"""

from stim_classes2 import *
from stream_methods2 import *

from pylab import *
import scipy.signal as signal

plt.cla
plt.clf
plt.close('all')

trigSamps = 1
CTRL_OFF = -1

silence5 = SigBlock('',{'tlen':5})
silence1 = SigBlock('',{'tLen':1})


off1 = PulseBlock({'hiT':1,'hiVal':CTRL_OFF}) 

LPcut = 200.0
filtP = {}#{'LPcutHz':LPcut}

varScale = 1#25e3/LPcut#125 #compensate for the 200Hz filtering


"""
STEP FR
do this from controller?
#
"""
"""
CTRL_step response
"""

'''
fold = "C_singlePulse/x/"

st = StimStream()

sb = PulseBlock({'hiT':1.5,'lowT':3.5,'hiVal':1})
offb = PulseBlock({'tStart':sb.hiT,'hiT':sb.lowT,'lowT':0,'hiVal':CTRL_OFF})
st.addBlock(sb)
st.addBlock(offb)
st = st.collapse()#st.compileBlocks()
st.HPLPArray(filtP)

st.assignFolder(fold)
st.printToTxt('C_singlePulseFR')
tr = st.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_singlePulseFR'})
'''

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#CTRL Charac
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

'''
fold1 = "C_wn_1x/x/"
fold2 = "C_wn_10x/x/"

noiseRange = 1

var = varScale*(noiseRange/6.0)**2
dc = .5
off5 = PulseBlock({'hiT':3.5,'hiVal':CTRL_OFF})

w = StimStream()
noiseB = WNBlock({'tLen':1.5,'mean':dc,'var':var,'lBound':0,'uBound':1})
w.addBlock(noiseB)
w.appendBlock(off5)
w.compileBlocks()
print 
print
w10 = repTrain(w,w.tLen,10)

#collapses the two blocks together
wc = w.collapse()
wc.HPLPArray(filtP)

w.assignFolder(fold1)
wc.assignFolder(fold1)
wc.printToTxt('C_wnDC(.5)_std(0p0277)_1rep')
tr = w.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_wnDC(.5)_std(0p0277)_1rep'})

w.plot()

print 'start replicating'
print ''



w10.compileBlocks()
w10.HPLPArray(filtP)
w10.assignFolder(fold2)
w10.printToTxt('C_wnDC(.5)_std(0p0277)_10rep')
tr = w10.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_wnDC(.5)_std(0p0277)_10rep'})
w10.plot()
'''

'''
#sines
'''
fold1 = 'C_sine1/x/'
fold2 = 'C_sine10/x/'
fold3 = 'C_sine20/x/'
fold4 = 'C_sine30/x/'
fold5 = 'C_sineSweep/x/'
off5 = PulseBlock({'hiT':3.5,'lowT':0,'hiVal':CTRL_OFF})
ampRange = 1.0
dc = .5

freqs = np.linspace(0,45,num=10)
freqs[0] = 1

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':1.5,'dcOff':dc,'freqHz':1,'amp':ampRange/2,'phase0':0}))
s.appendBlock(off5)
s.HPLPArray(filtP)
s.assignFolder(fold1)
s.printToTxt('C_sine_DC(.5)_amp(.5)_1Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(.5)_amp(.5)_1Hz'})

s.blocks[0].modTheseParams({'freqHz':10})
s.HPLPArray(filtP)
s.assignFolder(fold2)
s.printToTxt('C_sine_DC(.5)_amp(.5)_10Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(.5)_amp(.5)_10Hz'})


s.blocks[0].modTheseParams({'freqHz':20})
s.HPLPArray(filtP)
s.assignFolder(fold3)
s.printToTxt('C_sine_DC(.5)_amp(.5)_20Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(.5)_amp(.5)_20Hz'})

s.blocks[0].modTheseParams({'freqHz':30})
s.HPLPArray(filtP)
s.assignFolder(fold4)
s.printToTxt('C_sine_DC(.5)_amp(.5)_30Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(.5)_amp(.5)_30Hz'})
s.plot()

#sweep
sweep = StimStream()

for f in freqs:
    sweep.appendBlock(SineBlock({'tLen':1.5,'dcOff':dc,'freqHz':f,'amp':ampRange/2,'phase0':0}))
    sweep.appendBlock(off5)
sweep = sweep.collapse()
sweep.HPLPArray(filtP)
sweep.assignFolder(fold5)
sweep.printToTxt('C_sineSweep_DC(.5)_amp(.5)_1to45Hz_10Steps')
sweep.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sineSweep_DC(.5)_amp(.5)_1to45Hz_10Steps'})

