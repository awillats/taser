# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 15:32:03 2015

@author: Adam
"""

from stim_classes import *
from stream_methods import *

from pylab import *
import scipy.signal as signal

plt.cla
plt.clf
plt.close('all')

trigSamps = 1
CTRL_OFF = -1

silence5 = SigBlock('',0,5)
silence1 = SigBlock('',0,1)
off1 = PulseBlock({'tStart':0,'hiT':1,'lowT':0,'hiVal':CTRL_OFF}) 


"""
#OL Light responses
"""
#Steps
f='O_randSteps/x/'
os.makedirs(f)

nLevels = 15
randSteps = genBlockTrain('Pulse',0,4,{'hiT':1.5,'lowT':2.5,'hiVal':1},nLevels)
#random
randAmps = np.random.rand(nLevels)
#or pseudo
pseudoAmps = np.array(range(nLevels))/float(nLevels)
np.random.shuffle(pseudoAmps)

randSteps.modParam({'hiVal':randAmps})
randSteps.compileBlocks()
randSteps.printToTxt(f+'O_randSteps_15')
tr = randSteps.genTriggers(trigSamps,{'doPrint':True,'fName':f+'ST_randSteps_15'})

#plt.plot(randSteps.getTArray(),tr)
#randSteps.plot()

f='O_pseudo/x/'
os.makedirs(f)
randSteps.modParam({'hiVal':randAmps})
randSteps.compileBlocks()
randSteps.printToTxt(f+'O_pseudoSteps_15')
tr = randSteps.genTriggers(trigSamps,{'doPrint':True,'fName':f+'ST_pseudoSteps_15'})
"""
"""

f1h='O_wn_dchi_1x/x/'
f1l='O_wn_dclo_1x/x/'
f5h='O_wn_dchi_5x/x/'
f5l='O_wn_dclo_5x/x/'
os.makedirs(f1h)
os.makedirs(f1l)
os.makedirs(f5h)
os.makedirs(f5l)

#WhiteNoise+dc
noiseRange = .3
var = (noiseRange/6.0)**2

dclo = noiseRange/2
dchi = .5-noiseRange/2


loW = StimStream()
loW.addBlock(WNBlock({'tStart':0,'tLen':3,'mean':dclo,'var':var}))
loW.appendBlock(silence5)

hiW = StimStream()
hiW.addBlock(WNBlock({'tStart':0,'tLen':3,'mean':dchi,'var':var}))
hiW.appendBlock(silence5)

#collapses the two blocks together
loWC = loW.collapse()
hiWC = hiW.collapse()

loWC.printToTxt(f1l+'O_wnDC(0p15)_std(0p05)_1rep')
hiWC.printToTxt(f1h+'O_wnDC(0p35)_std(0p05)_1rep')
tr = loW.genTriggers(trigSamps,{'doPrint':True,'fName':f1l+'ST_wnDC(0p15)_std(0p05)_1rep'})
tr = hiW.genTriggers(trigSamps,{'doPrint':True,'fName':f1h+'ST_wnDC(0p35)_std(0p05)_1rep'})

lo5 = replicateTrain(loW,loW.tLen,5)
hi5 = replicateTrain(hiW,hiW.tLen,5)

lo5.compileBlocks()
hi5.compileBlocks()
lo5.printToTxt(f5l+'O_wnDC(0p15)_std(0p05)_5rep')
hi5.printToTxt(f5h+'O_wnDC(0p35)_std(0p05)_5rep')
tr = lo5.genTriggers(trigSamps,{'doPrint':True,'fName':f5l+'ST_wnDC(0p15)_std(0p05)_5rep'})
tr = hi5.genTriggers(trigSamps,{'doPrint':True,'fName':f5h+'ST_wnDC(0p35)_std(0p05)_5rep'})
"""


"""
#sines

f1='O_sine1/x/'
f8='O_sine8/x/'
f16='O_sine16/x/'
f32='O_sine32/x/'
os.makedirs(f1)
os.makedirs(f8)
os.makedirs(f16)
os.makedirs(f32)
ampRange = .3
dc = .5-ampRange/2.0

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':1,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
sc=s.collapse()
sc.printToTxt(f1+'O_sine_DC(0p35)_1Hz')
sc.genTriggers(trigSamps,{'doPrint':True,'fName':f1+'ST_sine_DC(0p35)_1Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':8,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
s=s.collapse()
s.printToTxt(f8+'O_sine_DC(0p35)_8Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':f8+'ST_sine_DC(0p35)_8Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':16,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
s=s.collapse()
s.printToTxt(f16+'O_sine_DC(0p35)_16Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':f16+'ST_sine_DC(0p35)_16Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':32,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
s=s.collapse()
s.printToTxt(f32+'O_sine_DC(0p35)_32Hz')
tr=s.genTriggers(trigSamps,{'doPrint':True,'fName':f32+'ST_sine_DC(0p35)_32Hz'})
"""
"""
#ramp shapes
f='O_pulseShapes/x/'
os.makedirs(f)

rampT = 3e-3
rub = RampBlock({'tStart':0,'tLen':rampT,'dir':'UPRIGHT','amp':1})
rdb = RampBlock({'tStart':0,'tLen':rampT,'dir':'UPLEFT','amp':1})
pb = PulseBlock({'tStart':0,'hiT':rampT,'lowT':0,'hiVal':1})
gb = ProbeBlock({'tCenter':rampT,'tWide':rampT*2,'amp':1,'std':rampT/6.0})




tri = StimStream()
tri.addBlock(rub)
tri.appendBlock(rdb)
tri.appendBlock(SigBlock('',0,1-2*rub.tLen))
tri = tri.collapse()

rhb = StimStream()
rhb.addBlock(rub)
rhb.appendBlock(pb)
rhb.appendBlock(SigBlock('',0,1-2*rub.tLen))
rhb = rhb.collapse()

spb = StimStream()
spb.addBlock(pb)
spb.appendBlock(SigBlock('',0,1-pb.tLen))
spb = spb.collapse()

gpb = StimStream()
gpb.addBlock(gb)
gpb.appendBlock(SigBlock('',0,1-gb.tLen))
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
print 'precom'    
r.compileBlocks()
r.printToTxt(f+'O_tri_rh_p_g_3ms_4s')
tr = r.genTriggers(trigSamps,{'doPrint':True,'fName':f+'ST_tri_rh_p_g_3ms_4s'})
