#! /usr/bin/python
import LoadTraces
import des_block
import numpy as np
import matplotlib.pyplot as plt
from util import *

class DES_CPA():

    def __init__(self, acq):
        self.acq = acq
        self.distList = None

    def distinguisher(self, sboxNum, startTrace = 0, endTrace = 0, samplePoint = None, leakage='HW'):
        """Returns distinguisher values for specific sbox as a numpy float array
        If samplePoint is provided, then the distinguisher uses the correlation coefficient at this sample Point
        If no samplePoint is provided, the distinguisher uses the maximum, absolute cc from all samples in that traces
        """
        self.distList = np.zeros(64, dtype = np.float)
        if endTrace == 0:
             endTrace = self.acq.numTraces

        for i in range(0, 64):

            if (samplePoint == None and self.acq.numTraces>1):
                cc = self.corrCoefficient(i, sboxNum, startTrace, endTrace, leakage)
                index = np.argmax(np.absolute(cc))     # non-fixed point
                maxCorr = cc[index]
            else:
                if (self.acq.numTraces == 1):
                    samplePoint = 1
                ccSp = self.corrCoefficientSinglePoint(i, sboxNum, samplePoint, startTrace, endTrace, leakage)
                maxCorr = ccSp         # fixed point

            self.distList[i] = maxCorr
            #print self.distList[0x1f]
            #print "for key %x, DoM is %f" % (i, maxAbs)
        return self.distList


    def corrCoefficient(self, subkey, sboxNum, startTrace = 0, endTrace = 0, leakage='HW'):
        """
        Calculate CPA correlation co-efficient:
        a. Calculate sbox output for the given subkey
        b. Calculate hypothesis (which is hamming weight or hamming distance depending on power model)
        c. Calculate pearson's correlation co-efficient of given subkey

        Input:
        subkey (6bit key which is input to sbox)
        sboxNum (1...8)

        Return the correlation co-efficient of the given subkey
        """
        # Based on ChipWhisperer CPASimpleLoop.py

        sumnum = np.zeros(self.acq.numSamples)
        sumden1 = np.zeros(self.acq.numSamples)
        sumden2 = np.zeros(self.acq.numSamples)

        if (endTrace == 0):
            endTrace = self.acq.numTraces

        # Generate hypotheticals
        hyp = np.zeros(endTrace - startTrace)
        for i in range(endTrace - startTrace):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxIn = ptBlock ^ subkey
            sboxOut = des_block.sbox(sboxNum, sboxIn)
            if leakage =='HW':
                hyp[i] = HW(sboxOut)
            elif leakage =='HDRound':
                hyp[i] = HD(sboxOut,des_block.ptRoundIn(ptStr,sboxNum))
            elif leakage =='HD':
                hyp[i]  = HD(sboxOut, sboxIn)
            else:
                print "undefined leakage model"
                return 0


        meanh = np.mean(hyp, dtype = np.float)
        meant = np.mean(self.acq.traces[startTrace:endTrace], axis=0, dtype=np.float)

        for tnum in range(endTrace - startTrace):
            hdiff = (hyp[tnum] - meanh)
            tdiff = self.acq.traces[tnum+startTrace,:] - meant

            sumnum += hdiff * tdiff
            sumden1 += hdiff * hdiff
            sumden2 += tdiff * tdiff
        #print sumden1
        #print sumden2
        #print np.sqrt(sumden1 * sumden2)
        corr = sumnum / np.sqrt (sumden1 * sumden2)
        return corr


    def corrCoefficientSinglePoint(self, subkey, sboxNum, samplePoint, startTrace = 0, endTrace = 0, leakage ='HW'):
        """
        Calculate CPA correlation co-efficient:
        a. Calculate sbox out given subkey
        b. Calculate hypothesis (which is hamming weight or hamming distance depending on power model)
        c. Calculate pearson's correlation co-efficient of given subkey

        Input:
        subkey (6bit key which is input to sbox)
        sboxNum (1...8)

        Return the correlation co-efficient of the given subkey
        """
        # Based on ChipWhisperer CPASimpleLoop.py

        sumnum = 0.0
        sumden1 = 0.0
        sumden2 = 0.0

        if (endTrace == 0):
            endTrace = self.acq.numTraces

        # Generate hypotheticals
        hyp = np.zeros(endTrace - startTrace)
        for i in range(endTrace - startTrace):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxIn = ptBlock ^ subkey
            sboxOut = des_block.sbox(sboxNum, sboxIn)
            if leakage =='HW':
                hyp[i] = HW(sboxOut)   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HDRound':
                hyp[i] = HD(sboxOut,des_block.ptRoundIn(ptStr,sboxNum))   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HD':
                hyp[i]  = HD(sboxOut, sboxIn)
            else:
                print "undefined leakage model"
                return 0


        meanh = np.mean(hyp, dtype = np.float)
        meant = np.mean(self.acq.traces[startTrace:endTrace, samplePoint], dtype=np.float)

        for tnum in range(endTrace - startTrace):
            hdiff = (hyp[tnum] - meanh)
            tdiff = self.acq.traces[tnum+startTrace,:][samplePoint] - meant

            sumnum += hdiff * tdiff
            sumden1 += hdiff * hdiff
            sumden2 += tdiff * tdiff
        corr = sumnum / ((sumden1 * sumden2) ** 0.5)
        return corr




    def corrCoefficientSinglePointRivian(self, subkey, sboxNum, samplePoint, startTrace = 0, endTrace = 0, leakage ='HW'):
        """
        Calculate CPA correlation co-efficient:
        a. Calculate sbox out given subkey
        b. Calculate hypothesis (which is hamming weight or hamming distance depending on power model)
        c. Calculate pearson's correlation co-efficient of given subkey

        Input:
        subkey (6bit key which is input to sbox)
        sboxNum (1...8)

        Return the correlation co-efficient of the given subkey
        """
        # Based on ChipWhisperer CPASimpleLoop.py

        if (endTrace == 0):
            endTrace = self.acq.numTraces
        coeff = []
        # Generate hypotheticals
        hyp = np.zeros(endTrace - startTrace)
        for i in range(endTrace - startTrace):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxIn = ptBlock ^ subkey
            sboxOut = des_block.sbox(sboxNum, sboxIn)
            if leakage =='HW':
                hyp[i] = HW(sboxOut)   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HDRound':
                hyp[i] = HD(sboxOut,des_block.ptRoundIn(ptStr,sboxNum))   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HD':
                hyp[i]  = HD(sboxOut, sboxIn)
            else:
                print "undefined leakage model"
                return 0
            coeff.append(hyp[i] * self.acq.traces[i + startTrace])
        corr = np.mean(coeff)
        return corr




    def ptNumpyArray2String(self, inputtext):
        ptStr = ""
        for i in range(self.acq.blockSize):
            ptStr += hex2str(inputtext[i])
        return ptStr

if __name__ == "__main__":
    acq = LoadTraces.LoadTraces('traces/1001.trs', 100000)

    cpa = DES_CPA(acq)
    sboxNum = 1
    corr = cpa.corrCoefficient(0x1f,sboxNum, leakage= 'HDRound')
    plt.figure()
    plt.plot(corr, label='correct key')
    plt.draw()
    plt.show()

