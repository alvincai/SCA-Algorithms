#! /usr/bin/python

import LoadTraces
import des_block
import numpy as np
import matplotlib.pyplot as plt
from util import *

class DES_DPA():

    def __init__(self, acq):
        self.acq = acq
        self.distList = None


    def distinguisher(self, sboxNum, startTrace = 0, endTrace = 0, samplePoint = None):
        """Returns distinguisher values for specific sbox as a numpy float array
        If samplePoint is provided, then the distinguisher uses the DoM at this sample Point
        If no samplePoint is provided, the distinguisher uses the maximum, absolute DoM from all samples in the traces
        """

        self.distList = np.zeros(64, dtype = np.float)
        if endTrace == 0:
            endTrace = self.acq.numTraces

        for i in range(0, 64):
            if (samplePoint == None):
                dom = self.dom(i, sboxNum, startTrace = startTrace, endTrace = endTrace)
                maxAbs = np.amax(np.absolute(dom))     # non-fixed point
            else:
                domSP = self.domSinglePoint(i, sboxNum, samplePoint, startTrace = startTrace, endTrace = endTrace)
                maxAbs = np.absolute(domSP)         # fixed point

            self.distList[i] = maxAbs
            #print "for key %x, DoM is %f" % (i, maxAbs)
        return self.distList


    def dom(self, subkey, sboxNum, startTrace=0, endTrace=0):
        """
        Calculate DPA DOM:
        a. Calculate LSB of sbox out given subkey
        b. average all the traces with LSB= one. Repeat for LSB= zero
        c. Calculate the difference of means

        Input:
        subkey (6bit key which is input to sbox)
        sboxNum (1...8)

        Return a trace of the difference of means
        """
        traceLSB0 = np.zeros(self.acq.numSamples, dtype = np.float64)
        traceLSB1 = np.zeros(self.acq.numSamples, dtype = np.float64)
        diff = np.zeros(self.acq.numSamples, dtype = np.float64)
        numLSB0 = 0
        numLSB1 = 0
        if (endTrace == 0):
            endTrace = self.acq.numTraces


        for i in range(startTrace, endTrace):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i])
            #print ("ptStr is " + ptStr)
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxOut = des_block.sbox(sboxNum, ptBlock ^ subkey)

            if ((sboxOut & 1) == 1):
                # Sum traces (traceLSB1) with LSB = 1
                numLSB1 += 1
                for j in range(self.acq.numSamples):
                    traceLSB1[j] += self.acq.traces[i][j]
            else:
                # Sum traces (traceLSB0) with LSB = 0
                numLSB0 += 1
                for j in range(self.acq.numSamples):
                    traceLSB0[j] += self.acq.traces[i][j]
        for j in range(self.acq.numSamples):
            traceLSB0[j] = traceLSB0[j] / numLSB0
            traceLSB1[j] = traceLSB1[j] / numLSB1
            diff[j] = traceLSB1[j] - traceLSB0[j]

        #print "numLSB0 = %d, numLSB1 = %d" % (numLSB0, numLSB1)

        return diff

    def domSinglePoint(self, subkey, sboxNum, samplePoint, startTrace=0, endTrace=0):
        """
        Calculate DPA DOM for a single point:
        a. Calculate LSB of sbox out given subkey
        b. average all the traces with LSB= one. Repeat for LSB= zero
        c. Calculate the difference of means

        Input:
        subkey (6bit key which is input to sbox)
        sboxNum (1...8)

        Return a trace of the difference of means
        """
        traceLSB0 = 0.0
        traceLSB1 = 0.0
        diff = 0.0
        numLSB0 = 0
        numLSB1 = 0
        if (endTrace == 0):
            endTrace = self.acq.numTraces


        for i in range(startTrace, endTrace):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i])
            #print ("ptStr is " + ptStr)
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxOut = des_block.sbox(sboxNum, ptBlock ^ subkey)

            if ((sboxOut & 1) == 1):
                # Sum traces (traceLSB1) with LSB = 1
                numLSB1 += 1.0
                traceLSB1 += self.acq.traces[i][samplePoint]
            else:
                # Sum traces (traceLSB0) with LSB = 0
                numLSB0 += 1.0
                traceLSB0 += self.acq.traces[i][samplePoint]
        traceLSB0 = traceLSB0 / numLSB0
        traceLSB1 = traceLSB1 / numLSB1
        diff = traceLSB1 - traceLSB0
        #print "traceLSB0 = %f, traceLSB1 = %f" %(traceLSB0, traceLSB1)
        #print "numLSB0 = %f, numLSB1 = %f" % (numLSB0, numLSB1)

        return diff

    def ptNumpyArray2String(self, inputtext):
        ptStr = ""
        for i in range(self.acq.blockSize):
            ptStr += hex2str(inputtext[i])
        return ptStr

if __name__ == "__main__":
    acq = LoadTraces.LoadTraces('traces/0002.trs')
    dpa = DES_DPA(acq)
    dom = dpa.dom(0x1f, 1,0, 1000)
    plt.plot(dom)
    plt.show()




