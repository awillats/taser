# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 22:45:04 2015

@author: Adam
"""

from operator import attrgetter
import numpy as np
from stim_classes import *
import copy

class StimStream(SigBlock):
    """
    .blocks
    .tStarts
    """
    def __init__(self):
        #print 'nope'
        SigBlock.__init__(self,'Stream',0,0,{})
        self.blocks = []
        self.tStarts = []
        
    def addBlock(self,newBlock):
        self.blocks.append(copy.deepcopy(newBlock))
        
    def appendBlock(self,newBlock):
        self.compileBlocks()
        
        nb = copy.deepcopy(newBlock)
        p = nb.params
        p['tStart'] = self.tEnd
        
        if nb.sigType == '' or nb.sigType == 'Static':
            nb.tStart = p['tStart']
        else:
            nb.setFromParams(p)
        
        self.blocks.append(nb)    
    def compileBlocks(self,overrideParams={'Nope':0}):
        #find span
        #http://stackoverflow.com/questions/18005172/get-the-object-with-the-max-attributes-value-in-a-list-of-objects
    
        #set the whole stream to be at the largest freq
        self.SAMPLE_FREQ = max(b.SAMPLE_FREQ for b in self.blocks)    
    
        firstB = min(self.blocks, key=attrgetter('tStart'))        
        lastB = max(self.blocks, key=lambda x: (x.params['tStart']+x.params['tLen']))
        
        self.tStart = firstB.tStart
        self.tEnd = lastB.tStart+lastB.tLen
        
        if not overrideParams.has_key('Nope'):
            if overrideParams.has_key('endT') and overrideParams['endT']>self.tEnd:
               self.tEnd = overrideParams['endT']
               print 'end time has been increased to :',self.tEnd
            if overrideParams.has_key('startT') and overrideParams['startT']<self.tStart:
               self.tStart = overrideParams['startT']
               print 'start time has been changed to :',self.tStart
               
        self.tLen = self.tEnd-self.tStart
        self.genArray() 
        
        for b in self.blocks:      
            if b.SAMPLE_FREQ != self.SAMPLE_FREQ:
                #golden
                print 'one of these blocks needs interpolating'
                b.SAMPLE_FREQ = self.SAMPLE_FREQ         
            
            if b.sigType != 'Static':
                b.genArray()
            
            print b.sigType
            self.scaleShiftArraySection(b.tStart,b.tLen,1,b.valArray)
            self.tStarts.append(b.tStart)
            
        #condense it into compiled block??
        compBlock = SigBlock('Static',self.tStart,self.tLen,b.params)
        compBlock.valArray = self.valArray
        return compBlock
    def collapse(self):
        s = StimStream()
        c = self.compileBlocks()
        s.addBlock(c)
        s.compileBlocks()
        return s
        
    def modParam(self,newParams):
        """
        allows you to change parameters of a string of blocks
        """
        for k in newParams.keys():
            try:
                L = len(newParams[k]) #if this works it's an array
                if L != len(self.blocks):
                    print 'length mismatch between num plocks and num params',L,'!=',len(self.blocks)
                for i, b in enumerate(self.blocks):
                    b.params[k] = newParams[k][i]
                    b.setFromParams(b.params)
            except:
                for b in self.blocks:
                    b.params[k]=newParams[k]
                    #print b.params
                    b.setFromParams(b.params)
    def genTriggers(self, nSampHi=1, printParams={'doPrint':False}):
        TrigTrain = StimStream()
        
        tHi = float(nSampHi)/self.SAMPLE_FREQ        
        
        for b in self.blocks:
            #by default use tCenter for trigger point of probe block
            tTrig = b.tCenter if b.sigType == 'Probe' else b.tStart
            TrigTrain.addBlock(PulseBlock({'tStart':tTrig,'hiT':tHi,'lowT':0,'hiVal':1}))

        self.trigArray = TrigTrain.compileBlocks({'startT':self.tStart,'endT':self.tEnd}).valArray
        
        if printParams['doPrint'] == True:
            np.savetxt(printParams['fName']+'.txt',self.trigArray)
        return self.trigArray                
            
    def shiftStarts(self,tShift):
        #FIX ME
        for b in self.blocks:
            #print ">",b.params
            b.params['tStart'] = tShift + copy.deepcopy(b.params['tStart'])           
            
            if b.sigType == '' or b.sigType == 'Static':
                b.tStart += tShift
            else:
                if b.sigType == 'Probe':
                    b.params['tCenter'] += tShift
                b.setFromParams(b.params)
            
    def markStarts(self):
         print "'mark starts' said the parrot"
         print self.tStarts
         #raster of the start times        
    def plotEach(self):
        for b in self.blocks:
            plt.plot(b.getTArray(),b.valArray)



def genBlockTrain(blockType,t0,tLen,params,numReps):
    train = StimStream()

    tStarts = np.linspace(t0,t0+tLen*(numReps-1),num=numReps)    
    
    bParams = copy.deepcopy(params)
    
    if 'tLen' not in params:
        bParams['tLen'] = tLen
    
    for t in tStarts:
        bParams['tStart'] = t
        b = chooseBlock(blockType,bParams)
        train.addBlock(copy.deepcopy(b))
        
    return train


def genBlockTrainFromTrain(refTrain, newBlockType, startOffset, lenOffset, params):
    numReps = len(refTrain.tStarts)
    print numReps
    
    thisP = copy.deepcopy(params)
    
    trainOut = StimStream()    
    
    for b in refTrain.blocks:
        thisP['tStart'] = b.tStart + startOffset
        thisP['tLen'] = b.tLen + lenOffset
        b = chooseBlock(newBlockType,thisP)
        trainOut.addBlock(b)
    return trainOut
    
    
def trainFromBlock(block,tLen,numReps):
    train = StimStream()
    tInc = tLen
    bClone = copy.deepcopy(block)
    
    #change this garbage later
            
    for i in range(numReps):
        train.addBlock(bClone)
        #deep copies prevent references from being passed
        bClone = copy.deepcopy(bClone) 
        bClone.params['tStart'] += tInc
        
        if bClone.sigType == 'Static':
            bClone.tStart = bClone.params['tStart']
        else:
            bClone.setFromParams(bClone.params)
    #train.plotEach()
    return train
def replicateTrain(trainUnit,tLen,numReps):
    train = StimStream()
    tInc = tLen
    uClone = copy.deepcopy(trainUnit)
    
    #change this garbage later
            
    for i in range(numReps):
        for b in uClone.blocks:
            train.addBlock(b)
        
        uClone = copy.deepcopy(uClone) 
        #shift all starts by tInc
        uClone.shiftStarts(tInc)        
        
    return train       
"""
method ideas

test text files'

better setFromParams

Xbulk parameter modification
--scale shift version



?print to file
Xfilter
Xplace around
XrandomSteps
XampedPulseTrain


Xbound

"""