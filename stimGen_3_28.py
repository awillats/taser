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
STEP FR
do this from controller?
#
"""
"""
CTRL_step response
"""

st = StimStream()

sb = PulseBlock({'tStart':0,'hiT':2,'lowT':3,'hiVal':1})
offb = PulseBlock({'tStart':sb.hiT,'hiT':sb.lowT,'lowT':0,'hiVal':CTRL_OFF})
st.addBlock(sb)
st.addBlock(offb)
st = st.collapse()#st.compileBlocks()
st.printToTxt('C_singlePulseFR')
tr = st.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_singlePulseFR'})


"""
#OL Light responses
"""
"""
#Steps
nLevels = 15
randSteps = genBlockTrain('Pulse',0,4,{'hiT':1.5,'lowT':2.5,'hiVal':1},nLevels)
#random
randAmps = np.random.rand(nLevels)
#or pseudo
pseudoAmps = np.array(range(nLevels))/float(nLevels)
np.random.shuffle(pseudoAmps)

randSteps.modParam({'hiVal':randAmps})
randSteps.compileBlocks()
randSteps.printToTxt('O_randSteps_15')
tr = randSteps.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_randSteps_15'})

#plt.plot(randSteps.getTArray(),tr)
#randSteps.plot()


randSteps.modParam({'hiVal':randAmps})
randSteps.compileBlocks()
randSteps.printToTxt('O_pseudoSteps_15')
tr = randSteps.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_pseudoSteps_15'})
"""
"""
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

loWC.printToTxt('O_wnDC(0p15)_std(0p05)_1rep')
hiWC.printToTxt('O_wnDC(0p35)_std(0p05)_1rep')
tr = loW.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_wnDC(0p15)_std(0p05)_1rep'})
tr = hiW.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_wnDC(0p35)_std(0p05)_1rep'})

lo5 = replicateTrain(loW,loW.tLen,5)
hi5 = replicateTrain(hiW,hiW.tLen,5)

lo5.compileBlocks()
hi5.compileBlocks()
lo5.printToTxt('O_wnDC(0p15)_std(0p05)_5rep')
hi5.printToTxt('O_wnDC(0p35)_std(0p05)_5rep')
tr = lo5.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_wnDC(0p15)_std(0p05)_5rep'})
tr = hi5.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_wnDC(0p35)_std(0p05)_5rep'})
"""


"""
#sines
ampRange = .3
dc = .5-ampRange/2.0

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':1,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
sc=s.collapse()
sc.printToTxt('O_sine_DC(0p35)_1Hz')
sc.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(0p35)_1Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':8,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
s=s.collapse()
s.printToTxt('O_sine_DC(0p35)_8Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(0p35)_8Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':16,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
s=s.collapse()
s.printToTxt('O_sine_DC(0p35)_16Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(0p35)_16Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':32,'amp':ampRange/2,'phase0':0}))
s.appendBlock(silence5)
s=s.collapse()
s.printToTxt('O_sine_DC(0p35)_32Hz')
tr=s.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_sine_DC(0p35)_32Hz'})
"""

#ramp shapes
rampT = 3e-3
rub = RampBlock({'tStart':0,'tLen':rampT,'dir':'UPRIGHT','amp':1})
rdb = RampBlock({'tStart':0,'tLen':rampT,'dir':'UPLEFT','amp':1})
pb = PulseBlock({'tStart':0,'hiT':rampT,'lowT':0,'hiVal':1})
gb = ProbeBlock({'tCenter':rampT,'tWide':rampT*2,'amp':1,'std':rampT/6.0})




tri = StimStream()
tri.addBlock(rub)
tri.appendBlock(rdb)
tri.appendBlock(PulseBlock({'tStart':0,'hiT':1-2*rub.tLen,'lowT':0,'hiVal':CTRL_OFF}))
tri = tri.collapse()

rhb = StimStream()
rhb.addBlock(rub)
rhb.appendBlock(pb)
rhb.appendBlock(PulseBlock({'tStart':0,'hiT':1-2*rub.tLen,'lowT':0,'hiVal':CTRL_OFF}))
rhb = rhb.collapse()

spb = StimStream()
spb.addBlock(pb)
spb.appendBlock(PulseBlock({'tStart':0,'hiT':1-pb.tLen,'lowT':0,'hiVal':CTRL_OFF}))
spb = spb.collapse()

gpb = StimStream()
gpb.addBlock(gb)
gpb.appendBlock(PulseBlock({'tStart':0,'hiT':1-gb.tLen,'lowT':0,'hiVal':CTRL_OFF}))
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
r.printToTxt('O_tri_rh_p_g_3ms_4s')
tr = r.genTriggers(trigSamps,{'doPrint':True,'fName':'ST_tri_rh_p_g_3ms_4s'})

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#CTRL Charac
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""
"""
noiseRange = 4
var = (noiseRange/6.0)**2
dc = 5

off5 = PulseBlock({'tStart':0,'hiT':5,'lowT':0,'hiVal':CTRL_OFF})

w = StimStream()
w.addBlock(WNBlock({'tStart':0,'tLen':3,'mean':dc,'var':var}))
w.appendBlock(off5)

#collapses the two blocks together
wc = w.collapse()

wc.printToTxt('C_wnDC(5)_std(0p333)_1rep')
tr = w.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_wnDC(5)_std(0p333)_1rep'})

w5 = replicateTrain(w,w.tLen,5)

w5.compileBlocks()
w5.printToTxt('C_wnDC(5)_std(0p333)_5rep')
tr = w5.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_wnDC(5)_std(0p333)_5rep'})


#sines
ampRange = 4
dc = 5

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':1,'amp':ampRange/2,'phase0':0}))
s.appendBlock(off5)
s=s.collapse()
s.printToTxt('C_sine_DC(5)_amp(pm2)_1Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_1Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':8,'amp':ampRange/2,'phase0':0}))
s.appendBlock(off5)
s=s.collapse()
s.printToTxt('C_sine_DC(5)_amp(pm2)_8Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_8Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':16,'amp':ampRange/2,'phase0':0}))
s.appendBlock(off5)
s=s.collapse()
s.printToTxt('C_sine_DC(5)_amp(pm2)_16Hz')
s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_16Hz'})

s = StimStream()
s.addBlock(SineBlock({'tStart':0,'tLen':3,'dcOff':dc,'freqHz':32,'amp':ampRange/2,'phase0':0}))
s.appendBlock(off5)
s=s.collapse()
s.printToTxt('C_sine_DC(5)_amp(pm2)_32Hz')
tr=s.genTriggers(trigSamps,{'doPrint':True,'fName':'CT_sine_DC(5)_amp(pm2)_32Hz'})
