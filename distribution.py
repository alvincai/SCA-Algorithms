#! /usr/bin/python
import LoadTraces
import des_block
import numpy as np
import matplotlib.pyplot as plt
from util import *

class Distribution():

    def __init__(self, acq):
        self.acq = acq


    def sboxInputs(self, subkey, sboxNum, startTrace =0, endTrace = 0):
        distribution = np.zeros(64)
        if (endTrace == 0):
            endTrace = self.acq.numTraces

        # Generate hypotheticals
        for i in range(endTrace - startTrace):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxIn = ptBlock ^ subkey
            distribution[sboxIn] +=1

        print distribution


    def ptNumpyArray2String(self, inputtext):
        ptStr = ""
        for i in range(self.acq.blockSize):
            ptStr += hex2str(inputtext[i])
        return ptStr

if __name__ == "__main__":
    acq = LoadTraces.LoadTraces('traces/0005.trs')
    d = Distribution(acq)
    print "sbox 2 of 0005"
    d.sboxInputs(0x31, 2)
