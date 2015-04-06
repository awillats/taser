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

silence5 = SigBlock('',{'tLen':5})
silence1 = SigBlock('',{'tLen':1})
off1 = PulseBlock({'tStart':0,'hiT':1,'lowT':0,'hiVal':CTRL_OFF}) 

filtP = {}

"""
#OL Light responses
"""

#Steps
f = 'O_singlePulse/x/'
b = StimStream()
b.addBlock(PulseBlock({'hiT':1.5,'lowT':3.5,'hiVal':1}))

b.plot()
b.assignFolder(f)
b.printToTxt('O_singlePulse_1p5_3p5')
b.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_singlePulse_1p5_3p5'})

''''
f='O_randSteps/x/'

nLevels = 15
randSteps = genBlockTrain('Pulse',0,4,{'hiT':1.5,'lowT':2.5,'hiVal':1},nLevels)
#random
randAmps = np.random.rand(nLevels)
#or pseudo
pseudoAmps = np.array(range(nLevels))/float(nLevels)
np.random.shuffle(pseudoAmps)

randSteps.modParam({'hiVal':randAmps})
randSteps.compileBlocks()
randSteps.assignFolder(f)
randSteps.printToTxt('O_randSteps_15')
tr = randSteps.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_randSteps_15'})

#plt.plot(randSteps.getTArray(),tr)
#randSteps.plot()

f='O_pseudo/x/'
randSteps.modParam({'hiVal':randAmps})
randSteps.compileBlocks()
randSteps.assignFolder(f)
randSteps.printToTxt('O_pseudoSteps_15')
tr = randSteps.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_pseudoSteps_15'})
'''
"""
WHITE
"""
'''
fold1 = "O_wn_1x/x/"
fold2 = "O_wn_5x/x/"

noiseRange = 1

var = varScale*(noiseRange/6.0)**2
dc = .5
off5 = PulseBlock({'hiT':3.5,'hiVal':0})

w = StimStream()
noiseB = WNBlock({'tLen':1.5,'mean':dc,'var':var,'lBound':0,'uBound':1})
w.addBlock(noiseB)
w.appendBlock(off5)
w.compileBlocks()
w10 = repTrain(w,w.tLen,5)

#collapses the two blocks together
wc = w.collapse()
wc.HPLPArray(filtP)

w.assignFolder(fold1)
wc.assignFolder(fold1)
wc.printToTxt('O_wnDC(.5)_std(0p0277)_5rep')
tr = w.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_wnDC(.5)_std(0p0277)_5rep'})

w.plot()

print 'start replicating'
print ''



w10.compileBlocks()
w10.HPLPArray(filtP)
w10.assignFolder(fold2)
w10.printToTxt('O_wnDC(.5)_std(0p0277)_10rep')
tr = w10.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_wnDC(.5)_std(0p0277)_10rep'})
w10.plot()
'''



"""
Sines
"""

'''
#sines
'''

fold1 = 'O_sine1/x/'
fold2 = 'O_sine10/x/'
fold3 = 'O_sine20/x/'
fold4 = 'O_sine30/x/'
fold5 = 'O_sineSweep6/x/'
off5 = PulseBlock({'hiT':3.5,'lowT':0,'hiVal':0})
ampRange = 1.0
dc = .5

#freqs = np.linspace(0,45,num=10)
freqs = np.array([1,5,10,20,30,40])
freqs[0] = 1
'''
s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':1.5,'dcOff':dc,'freqHz':1,'amp':ampRange/2,'phase0':0}))
s.appendBlock(off5)
s.HPLPArray(filtP)
s.assignFolder(fold1)
s.printToTxt('O_sine_DC(.5)_amp(.5)_1Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(.5)_amp(.5)_1Hz'})

s.blocks[0].modTheseParams({'freqHz':10})
s.HPLPArray(filtP)
s.assignFolder(fold2)
s.printToTxt('O_sine_DC(.5)_amp(.5)_10Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(.5)_amp(.5)_10Hz'})


s.blocks[0].modTheseParams({'freqHz':20})
s.HPLPArray(filtP)
s.assignFolder(fold3)
s.printToTxt('O_sine_DC(.5)_amp(.5)_20Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(.5)_amp(.5)_20Hz'})

s.blocks[0].modTheseParams({'freqHz':30})
s.HPLPArray(filtP)
s.assignFolder(fold4)
s.printToTxt('O_sine_DC(.5)_amp(.5)_30Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(.5)_amp(.5)_30Hz'})
s.plot()
'''
#sweep
sweep = StimStream()

for f in freqs:
    sweep.appendBlock(SineBlock({'tLen':1.5,'dcOff':dc,'freqHz':f,'amp':ampRange/2,'phase0':0}))
    sweep.appendBlock(off5)
sweep = sweep.collapse()
sweep.HPLPArray(filtP)
sweep.assignFolder(fold5)
sweep.printToTxt('O_sineSweep_DC(.5)_amp(.5))_1to40Hz_6Steps')
sweep.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sineSweep_DC(.5)_amp(.5))_1to40Hz_6Steps'})



"""
#ramp shapes
"""
'''
f='O_pulseShapes/x/'

rampT = 3e-3
rub = RampBlock({'tStart':0,'tLen':rampT,'dir':'UPRIGHT','amp':1})
rdb = RampBlock({'tStart':0,'tLen':rampT,'dir':'UPLEFT','amp':1})
pb = PulseBlock({'tStart':0,'hiT':rampT,'lowT':0,'hiVal':1})
gb = ProbeBlock({'tCenter':rampT,'tWide':rampT*2,'amp':1,'std':rampT/6.0})

tri = StimStream()
tri.addBlock(rub)
tri.appendBlock(rdb)
tri.appendBlock(SigBlock('',{'tLen':1-2*rub.tLen}))
tri = tri.collapse()

rhb = StimStream()
rhb.addBlock(rub)
rhb.appendBlock(pb)
rhb.appendBlock(SigBlock('',{'tLen':1-2*rub.tLen}))
rhb = rhb.collapse()

spb = StimStream()
spb.addBlock(pb)
spb.appendBlock(SigBlock('',{'tLen':1-pb.tLen}))
spb = spb.collapse()

gpb = StimStream()
gpb.addBlock(gb)
gpb.appendBlock(SigBlock('',{'tLen':1-gb.tLen}))
gpb = gpb.collapse()


r = StimStream()
    
rhb.shiftStarts(1)
spb.shiftStarts(2)
gpb.shiftStarts(3)

for b in tri.blocks:
    r.addBlock(b)
for b in rhb.blocks:
    r.addBlock(b)
for b in spb.blocks:    
    r.addBlock(b)
for b in gpb.blocks:
    r.addBlock(b)
    
r.compileBlocks()
r.assignFolder(f)
r.printToTxt('O_tri_rh_p_g_3ms_4s')
tr = r.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_tri_rh_p_g_3ms_4s'})
'''