#! /usr/bin/python
import numpy as np
import LoadTraces
import DES_DPA
import DES_CPA
import des_block
from util import *
import pickle
import Confusion

class RDM():
    """Calculates the relative distinguishing margin as described in A Fair Evaluation Framework for Comparing Side-Channel
    Distinguishers by Carolyn Whitnall and Elisabeth Oswald"""

    def __init__(self, distinguisher, correctSubKey):
        """input distinguisher values as numpy array, and the correct SubKey value"""
        self.distinguisher = distinguisher
        self.correctSubKey = correctSubKey

    def rdm(self):
        self.distinguisher = np.absolute(self.distinguisher)
        self.distinguisherCorrectSubKey = self.distinguisher[self.correctSubKey]
        stddev = np.std(self.distinguisher)
        self.distinguisher[self.correctSubKey] = 0           # remove distinguisher value of correct key
        maxNotCorrect = np.amax(np.absolute(self.distinguisher))    # Obtain the maximum distinguisher value, which is not the correct key

        #print "correct subkey is %f" % self.distinguisherCorrectSubKey
        #print "maximum distg valyue which is not the correct key is %f" % maxNotCorrect
        #print "std dev %f" % stddev

        return (self.distinguisherCorrectSubKey - maxNotCorrect)/ stddev

    def estimator(self, signal, noise,numTraces):
        cov, mean = self.conf.YCovMeanCPA2(self.signal, self.noise, numTraces, correctKey=0)



    #def estimator(self,signal, correctKey=0, leakage='HW', sboxNum=1):
        #Nk = 64
        #numPT = 2 ** 6
        #cov = np.zeros(64)
        #ev = 2.0    # E [V | kc.] = E[V | kg] = 2.0

        #for i in range(64):
            #tmp = []
            #for ptBlock in range(numPT):
                #sboxOutc = des_block.sbox(sboxNum, ptBlock ^ correctKey)
                #sboxOuti = des_block.sbox(sboxNum, ptBlock ^ i)
                #if leakage =='HW':
                    #vkc = HW(sboxOutc) # V | kc
                    #vki = HW(sboxOuti) # V | ki
                #elif leakage =='HD':
                    #vkc = HD(sboxOutc, correctKey^ptBlock) # V | kc
                    #vki = HD(sboxOuti, i^ptBlock) # V | ki


                #tmp.append((vkc - 2.0) * (vki - 2.0) * signal)
            #cov[i] = np.mean(tmp)

        #return cov

if __name__ == "__main__":
    acq = LoadTraces.LoadTraces('traces/0002R.trs')

    # dpa = DES_DPA.DES_DPA(acq)
    # distg = dpa.distinguisher(1)

    #cpa = DES_CPA.DES_CPA(acq)
    #distg = cpa.distinguisher(1)

    #print "Distinguisher values are: ",
    #print distg

    rdm = RDM([1,1], 0x0)
    corr = []
    for sk in range(64):
        corr.append(rdm.estimator(1, sk))
    #print "Relative Distinguishing Metric for sbox1 is: ",
    #print rdm.rdm()


    pickle.dump(corr, open('correlation.p', 'wb'))

