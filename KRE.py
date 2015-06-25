#! /usr/bin/python

import des_block
import numpy as np
from matplotlib.pyplot import *
from scipy.stats import mvn
import scipy.io
import Confusion
import DES_DPA
import DES_CPA
import LoadTraces
from numpy.linalg import inv    # Used to calculate the inverse of the covariance matrix
from scipy.linalg import sqrtm  # Used to calculate the sqrt of the covariance matrix
import pickle

class KRE():

    """Class for calculating the Key Rank Evolution (Also known as partial guessing entropy) """

    def __init__(self, acq, sboxNum, correctSubKey,leakage='HW'):
        self.acq = acq
        self.sboxNum = sboxNum
        self.correctSubKey = correctSubKey
        self.leakage = leakage
        self.dpa = DES_DPA.DES_DPA(acq)
        self.cpa = DES_CPA.DES_CPA(acq)



    def rank(self, numTraces, samplePoint = None, attack='cpa', runs=1000, leakage='HW'):
        """
        Returns Empirical success rate for DPA for the indicated Number of measurements (numTraces)
        If samplePoint is provided, only uses single samplePoint, else use entire trace to calculate distinguishers
        Performs DPA attack 1000 times and returns the number of rank0 success/1000
        """
        rank = 0.0
        if ((numTraces * runs) > self.acq.numTraces):
            print "Not enough traces to calculate Empirical Success Rate. Need %d traces, have only %d traces" %((numTraces *runs), self.acq.numTraces)
            return 0.0

        if (attack =='cpa'):
            distinguisher = self.cpa.distinguisher
        else:
            distinguisher = self.dpa.distinguisher

        for i in range(runs):
            start = i*numTraces
            distg = distinguisher(self.sboxNum, startTrace = start, endTrace =(start + numTraces), samplePoint = samplePoint, leakage=leakage)
            distg = np.abs(distg)
            distCorrect = distg[self.correctSubKey]

            distg[self.correctSubKey] = 0
            for i in range(64):
                if distg[i] > distCorrect:
                    rank += 1.0

        return rank/runs


if __name__ == "__main__":
    trsFile = 'traces/0001.trs'
    acq = LoadTraces.LoadTraces(trsFile, samplePoint = 3881)
    correctKey = "DEADBEEFCAFEBABE"
    sboxNum = 1
    leakage = 'HW'
    roundNum = 1    # Do not Change this. Currently fixed to 1.
    rk = des_block.roundkey(correctKey,roundNum)   # Correct Round Key of roundNum
    sk = des_block.subkey(rk,sboxNum)

    r = KRE(acq, sboxNum, sk, leakage = leakage)
    yKRE = []
    xKRE = []
    for i in range (15, 17, 2):
        yKRE.append(r.rank(i, samplePoint = 0, attack='cpa', runs=10, leakage=leakage))
        xKRE.append(i)

    print yKRE
    print xKRE

