# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 22:45:04 2015

@author: Adam
"""

from operator import attrgetter
import numpy as np
from stim_classes2 import *
import copy

class StimStream(SigBlock):
    """
    .blocks
    .tStarts
    """
    #DEF_P = {'tStart':0,'tLen':0}
    def __init__(self,p={}):
        #print 'nope'
        SigBlock.__init__(self,'Stream',p)
        self.blocks = []
        self.tStarts = []
        
    def addBlock(self,newBlock):
        self.blocks.append(copy.deepcopy(newBlock))
    def mergeStream(self,streamIn):
        for b in streamIn.blocks:
            self.addBlock(b)
    def appendBlock(self,newBlock):
        if not self.blocks:
            self.addBlock(newBlock)
            return
            
        self.compileBlocks()
        
        nb = copy.deepcopy(newBlock)
        p = nb.params
        p['tStart'] = self.tEnd
        
        #FIX should be able to remove if
        """"
        if nb.sigType == '' or nb.sigType == 'Static':
            nb.tStart = p['tStart']
        else:
            nb._updateFromP(p)
        nb.params = p 
        """
        nb.modTheseParams(p)
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
               #print 'end time has been increased to :',self.tEnd
            if overrideParams.has_key('startT') and overrideParams['startT']<self.tStart:
               self.tStart = overrideParams['startT']
               #print 'start time has been changed to :',self.tStart
               
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
        comP = copy.deepcopy(b.params)
        comP['tStart'] = self.tStart
        comP['tLen'] = self.tLen
        
        compBlock = SigBlock('Static',comP)
        compBlock.valArray = self.valArray
        return compBlock
    def collapse(self):
        s = StimStream()
        c = self.compileBlocks()
        s.addBlock(c)
        s.compileBlocks()
        s.sigType = 'Static'
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
                    b._updateFromP(b.params)
            except:
                for b in self.blocks:
                    b.params[k]=newParams[k]
                    #print b.params
                    b._updateFromP(b.params)
    def genTriggers(self, nSampHi=1, printParams={'doPrint':False}):
        print printParams
        print 1,4,5,6
        TrigTrain = StimStream()
        tHi = float(nSampHi)/self.SAMPLE_FREQ        
        
        for b in self.blocks:
            #by default use tCenter for trigger point of probe block
            tTrig = b.tCenter if b.sigType == 'Probe' else b.tStart
            TrigTrain.addBlock(PulseBlock({'tStart':tTrig,'hiT':tHi,'lowT':0,'hiVal':1}))
       
        self.trigArray = TrigTrain.compileBlocks({'startT':self.tStart,'endT':self.tEnd}).valArray
        
        if printParams['doPrint'] == True:
            if hasattr(self,'foldName'):
                printParams['fName'] = self.foldName+printParams['fName']
            np.savetxt(printParams['fName']+'.txt',self.trigArray)
            
        return self.trigArray
                
    def plotTriggers(self,nSampHi=1,printParams={'doPrint':False}):
        if not hasattr(self,'trigArray'):
            print 'triggers generated'
            genTriggers(nSampHi,printParams)
        plt.plot(self.getTArray(),self.trigArray)
        return self.trigArray
    def shiftStarts(self,tShift):
        for b in self.blocks:
            b.shiftStart(tShift)
    def markStarts(self):
         print "'mark starts' said the parrot"
         print self.tStarts
         #raster of the start times
    def plot(self):
        self.compileBlocks()
        SigBlock.plot(self)
        
    def plotEach(self):
        for b in self.blocks:
            b.plot()



def genBlockTrain(blockType,t0,tLen,params_,numReps):
    #could do this the same way as the others
    train = StimStream()
    tStarts = np.linspace(t0,t0+tLen*(numReps-1),num=numReps)    
    params = copy.deepcopy(params_)
    #if 'tLen' not in params:
     #   bParams['tLen'] = tLen
    
    for t in tStarts:
        params['tStart'] = t
        b = chooseBlock(blockType,params)
        train.addBlock(copy.deepcopy(b))
        
    return train
    
def repBlock(block,tLen,numReps):
    train = StimStream()
    tInc = tLen
    bClone = copy.deepcopy(block)
            
    for i in range(numReps):
        train.addBlock(bClone)
        #deep copies prevent references from being passed
        bClone = copy.deepcopy(bClone) 
        bClone.shiftStart(tInc)
        
    return train

def genTrainFromTrain(refTrain, newBlockType, startOffset, lenOffset, params):
    #numReps = len(refTrain.tStarts)
    #print numReps
    thisP = copy.deepcopy(params)
    trainOut = StimStream()    
    
    for br in refTrain.blocks:
        bo = chooseBlock(newBlockType,thisP)
        bo.modTheseParams({'tStart':br.tStart+startOffset,'tLen':br.tLen+lenOffset})
        trainOut.addBlock(bo)
    return trainOut
    
def repTrain(trainUnit,tLen,numReps):
    train = StimStream()
    
    tInc = tLen
    print 'tinc',tInc
    uClone = copy.deepcopy(trainUnit)
            
    for i in range(numReps):
        for b in uClone.blocks:
            train.addBlock(b)
        
        uClone = copy.deepcopy(uClone) 
        uClone.shiftStarts(tInc)   
    return train  
     
"""
method ideas
check for overlap when adding together
compile=plot?=print?

XplotTrigger
~better setFromParams
Xdefault folder to print to
XAPPEND BLOCK
Xtest text files'
?print to file
Xfilter
Xplace around
XrandomSteps
XampedPulseTrain

Xbulk parameter modification
--scale shift version

Xbound

"""