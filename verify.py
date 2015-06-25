#! /usr/bin/python
import LoadTraces
import des_block
import numpy as np
import matplotlib.pyplot as plt
from util import *
import distribution
from Crypto.Cipher import DES

class Verify():

    def __init__(self, acq):
        self.acq = acq
        self.distList = None
        self.key = des_block.des_block("deadbeefcafebabe", 64)
        self.des = DES.new('\xde\xad\xbe\xef\xca\xfe\xba\xbe', DES.MODE_ECB)

    def verify2(self):
        """
        Verify that ciphertext = Encrypt(plaintext, key) and delete those which fail this test (due to comms error)
        Uses the python DES Cipher block so is faster than verify()
        """
        count = 0
        deleteMe = []
        # Generate hypotheticals
        for i in range(self.acq.numTraces):
            ptStr = ''
            ctStr = ''
            for c in self.acq.inputtext[i]:
                ptStr += chr(c)
            for c in self.acq.outputtext[i]:
                ctStr += chr(c)
            ctStr1 = self.des.encrypt(ptStr)
            if (ctStr1 != ctStr):
                deleteMe.append(i)
                count += 1
                #print count
        print "Deleted %d traces" % count

        for i in reversed(deleteMe):
            self.acq.traces = np.delete(self.acq.traces, i, axis=0)
            self.acq.inputtext = np.delete(self.acq.inputtext, i, axis = 0)
            self.acq.outputtext = np.delete(self.acq.outputtext, i , axis =0)
        self.acq.numTraces -=  count


    def verify(self):
        """
        Verify that ciphertext = Encrypt(plaintext, key) and delete those which fail this test (due to comms error)
        """
        count = 0
        # Generate hypotheticals
        for i in range(self.acq.numTraces):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i])
            ctStr1 = (des_block.des_block(ptStr,64)).encipher(self.key)
            ctStr2 = self.ptNumpyArray2String(self.acq.outputtext[i])
            if (ctStr1 != ctStr2):
                acq.traces = np.delete(acq.traces, i, axis=0)
                acq.inputtext = np.delete(acq.inputtext, i, axis = 0)
                acq.outputtext = np.delete(acq.outputtext, i , axis =0)
                count += 1
        print "Deleted %d traces" % count



    def ptNumpyArray2String(self, inputtext):
        ptStr = ""
        for i in range(self.acq.blockSize):
            ptStr += hex2str(inputtext[i])
        return ptStr

if __name__ == "__main__":
    #acq = LoadTraces.LoadTraces('cw_tc6+trim.trs')
    acq = LoadTraces.LoadTraces('traces/0005.trs', 100000)
    v = Verify(acq)
    v.verify2()
    d = distribution.Distribution(acq)
    print "sbox 1 of 0005"
    d.sboxInputs(0x1f, 1)

