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

sb = PulseBlock({'hiT':2,'lowT':3,'hiVal':1})
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


fold1 = "C_wn_1x/x/"
fold2 = "C_wn_5x/x/"

noiseRange = 4

var = varScale*(noiseRange/6.0)**2
dc = 5
print 'off53'
print ''
print ''
off53 = PulseBlock({'hiT':5,'hiVal':CTRL_OFF})

w = StimStream()
w.addBlock(WNBlock({'tLen':3,'mean':dc,'var':var}))
w.appendBlock(off53)

#collapses the two blocks together
wc = w.collapse()
wc.HPLPArray(filtP)

w.assignFolder(fold1)
wc.assignFolder(fold1)
wc.printToTxt('C_wnDC(5)_std(0p333)_1rep')
tr = w.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_wnDC(5)_std(0p333)_1rep'})

print 'start replicating'
print ''

#w5 = trainFromBlock(wc,wc.tLen,5)
w5 = repTrain(w,w.tLen,5)

w5.compileBlocks()
w5.HPLPArray(filtP)
w5.assignFolder(fold2)
w5.printToTxt('C_wnDC(5)_std(0p333)_5rep')
tr = w5.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_wnDC(5)_std(0p333)_5rep'})
w5.plot()

#sines

fold1 = 'C_sine1/x/'
fold2 = 'C_sine8/x/'
fold3 = 'C_sine16/x/'
fold4 = 'C_sine32/x/'
off53 = PulseBlock({'tStart':3,'hiT':5,'lowT':0,'hiVal':CTRL_OFF})
ampRange = 4
dc = 5

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':1,'amp':ampRange/2,'phase0':0}))
s.addBlock(off53)
s=s.collapse()
s.HPLPArray(filtP)
s.assignFolder(fold1)
s.printToTxt('C_sine_DC(5)_amp(pm2)_1Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_1Hz'})


s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':8,'amp':ampRange/2,'phase0':0}))
s.addBlock(off53)
s=s.collapse()
s.HPLPArray(filtP)
s.assignFolder(fold2)
s.printToTxt('C_sine_DC(5)_amp(pm2)_8Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_8Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':16,'amp':ampRange/2,'phase0':0}))
s.addBlock(off53)
s=s.collapse()
s.HPLPArray(filtP)
s.assignFolder(fold3)
s.printToTxt('C_sine_DC(5)_amp(pm2)_16Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_16Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':32,'amp':ampRange/2,'phase0':0}))
s.addBlock(off53)
s=s.collapse()
s.HPLPArray(filtP)
s.assignFolder(fold4)
s.printToTxt('C_sine_DC(5)_amp(pm2)_32Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_32Hz'})

