# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 11:36:54 2015

@author: Adam
"""

import stim_classes2 as stc
import stream_methods2 as strm
import matplotlib.pyplot as plt
plt.cla
plt.close('all')
trigSamps = 1
"""
s = stc.SigBlock('',{'tStart':0,'tLen':5})

w = stc.WNBlock({'tStart':0,'tLen':2})
w = stc.WNBlock()
p=stc.PulseBlock({'hiT':5,'hiVal':-1})
#p.plot()

pr = stc.ProbeBlock()
p5 = strm.repBlock(pr,pr.tLen*5,6)
#p5.plotEach()
"""
q5 = strm.genBlockTrain('Ramp',0,3,{'dir':'asd'},7)
c = strm.genTrainFromTrain(q5,'Pulse',-.1,.2,{'hiT':.5,'lowT':0})

#c.mergeStream(q5)
#c.compileBlocks()
c.plotEach()

plt.figure()
ct = strm.repTrain(c,c.tLen+50,4)

#ct.plot first works
ct.plot()
#ct.plotEach()