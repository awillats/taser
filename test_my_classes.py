# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 20:26:47 2015

@author: Adam
"""
from stim_classes import *
from stream_methods import *

from pylab import *
import scipy.signal as signal

plt.cla
plt.clf
plt.close('all')

x = WNBlock({'tStart':0,'tLen':1,'mean':0,'var':1})
"""
#y = WNBlock(10,5+25e-6,0,2)
y = WNBlock({'tStart':10,'tLen':5,'minVal':0,'maxVal':1})
x.setRandSeed(123)

a = ProbeBlock({'tCenter':0,'tWide':10,'amp':2.5,'std':.5})
b = WNBlock({'tStart':30,'tLen':.1,'minVal':-1,'maxVal':-2})

ss = StimStream()
ss.addBlock(x)
ss.addBlock(y)
ss.addBlock(a)
ss.addBlock(b)
ss.compileBlocks()
"""
#plot(ss.getTArray(),ss.valArray)

"""
segment = StimStream()
segment.addBlock(WNBlock({'tStart':3,'tLen':1,'minVal':-1,'maxVal':1}))
segment.appendBlock(SigBlock('',0,5,{'tStart':0,'x':1}))
segment.appendBlock(RampBlock({'tStart':0,'tLen':5,'amp':.5,'dir':'UPRIGHT'}))
segBlock = segment.compileBlocks()

#segBlock.plot()
"""

#plt.figure()
#T = genBlockTrain('Pulse',0,100,{'hiT':10,'lowT':95,'hiVal':1.3},6)
#T = genBlockTrain('WN',0,100,{'tLen':5,'minVal':-1,'maxVal':1},6)
#T = trainFromBlock(segBlock,segBlock.tLen+1,6)
amps = np.random.rand(6)
amps = np.array([1,2,3,4,5,6])

#T = trainFromBlock(x,x.tLen*2,6)
#T.modParam({'maxVal':amps,'minVal':-amps})
#Tc = T.compileBlocks()
#Tc.boundArray(-.5,.5)
#Tc.plot()


#Q = genBlockTrainFromTrain(T,'Pulse',-.1,.3,{'hiVal':.6})

sb = StimStream()
b1 = SineBlock({'tStart':0,'tLen':2,'amp':1,'freqHz':1,'phase0':np.random.rand()*2*pi,'dcOff':0})
b2 = SineBlock({'tStart':0,'tLen':2,'amp':1,'freqHz':8,'phase0':np.random.rand()*2*pi,'dcOff':0})
b3 = SineBlock({'tStart':0,'tLen':2,'amp':1,'freqHz':16,'phase0':np.random.rand()*2*pi,'dcOff':0})
b4 = SineBlock({'tStart':0,'tLen':2,'amp':1,'freqHz':500,'phase0':np.random.rand()*2*pi,'dcOff':0})
b5 = SineBlock({'tStart':0,'tLen':2,'amp':1,'freqHz':1000,'phase0':np.random.rand()*2*pi,'dcOff':0})

bp = ProbeBlock({'tCenter':0.05,'tWide':2,'amp':10,'std':0.005})

sb.addBlock(b1)

#sb.addBlock(b2)
#sb.addBlock(b3)
sb.addBlock(b4)
sb.addBlock(b5)
#sb.addBlock(bp)
sbc = sb.compileBlocks()
#sbc.plot()
"""
Q = trainFromBlock(sbc,sbc.tLen*1.3,3)
Qc = Q.compileBlocks()
"""
print 'Q2'
Q2 = replicateTrain(sb,sb.tLen*1.1,3)
Q2c = Q2.compileBlocks()

print 'xt'
x = Q2.valArray#Tc.valArray
t = Q2.genTriggers(200)

plt.figure()
plt.subplot(2,1,1)
plt.plot(x)

f = Q2.HPLPArray({'n':2000,'LPcutHz':200})
plt.plot(f,'r')
plt.subplot(2,1,2)
plt.plot(t)



