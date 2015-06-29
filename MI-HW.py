#! /usr/bin/python

import LoadTraces
import des_block
import numpy as np
import matplotlib.pyplot as plt
from util import *
import scipy.stats
from collections import Counter
import pickle
from operator import mul
from fractions import Fraction

def nCk(n,k):
    return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1))

def entropyX():
    H_x = 0.0
    total = 0.0
    for i in range(5):
        total += nCk(4,i)
    for i in range(5):
        Pr_x = float(nCk(4,i)/total)
        H_x -= Pr_x * np.log2(Pr_x)
    return H_x

class DES_PI():
    """ Calculate the percieved information for DES """

    NK = 5                  # Number of DES Keys for an sbox
    H_x = 2.03063906223     # Entropy of sbox output in HD Round model
    Pr_x = [0.0625, 0.25, 0.375, 0.25, 0.0625]         # Probability of HD =0, HD = 1 ...

    def __init__(self, acq):
        self.acq = acq

        # For each key, the template consists of a mean, and co-variance matrix of noise.
        # Because we use a single point, n = 1, instead of a covariance matrix, we have a variance matrix. We save only the std dev.
        self.mean = np.zeros(DES_PI.NK)
        self.sd = np.zeros(DES_PI.NK)


    def pi(self, correctSubKey, sboxNum, samplePoint, startTrace = 0, endTrace = 0):
        """
        Returns perceived information value. Need to have first generated the template using buildTemplate function.

        # Inputs:
        correctSubKey - the correct subkey (which is 0 to 63 for DES)
        sboxNum - the sbox number we are targeting (values 1 to 8)
        samplePoint - the single sample point to calculate the PI for
        startTrace, endTrace - SNR is calculated from Trace startTrace to endTrace. This value is optional. If left out, it will calculate SNR for entire trace set
        """

        Nm = endTrace - startTrace
        pi = DES_PI.H_x

        leakage=[]
        for i in range(startTrace, endTrace):
            leakage.append(self.acq.traces[i][samplePoint])
        leakage = set(leakage)


        for i in range(int(min(leakage)-self.sd[0]*5), int(max(leakage)+self.sd[0]*5)):
            norm = 0.0
            for sboxCandidate in range(DES_PI.NK):
                # sum of PR[L|V] for all V
                norm += scipy.stats.norm(self.mean[sboxCandidate], self.sd[sboxCandidate]).pdf(i)  * DES_PI.Pr_x[sboxCandidate]

            for sboxIn in range(DES_PI.NK):
                PrLX = scipy.stats.norm(self.mean[sboxIn], self.sd[sboxIn]).pdf(i) #Pr[Leakage| V ]
                PrXL = PrLX  * DES_PI.Pr_x[sboxIn] / norm

                if(PrXL!=0.0):
                    if np.isnan(PrXL):
                        pass
                    else:
                        pi += DES_PI.Pr_x[sboxIn] * PrLX * np.log2(PrXL)
        return pi

    def check(self):
        for i in self.mean:
            if np.isnan(i):
                print "Error, not enough values collected"

    def buildTemplate(self, subkey, sboxNum, samplePoint, startTrace = 0, endTrace = 0, excludeStartTrace =None, excludeEndTrace=None):
        """
        Calculate template for each sbox input

        b. Calculate variance of the noise

        # Inputs:
        subkey - the correct subkey
        sboxNum - the sbox number we are targeting (values 1 to 8)
        samplePoint - the single sample point to calculate the PI for
        startTrace, endTrace - SNR is calculated from Trace startTrace to endTrace. This value is optional. If left out, it will calculate SNR for entire trace set
        excludestartTrace, excludedEndTrace - Exclude the traces between excludeStartTrace to excludeEndTrace. This interval should be between startTrace and endTrace. It is used for cross validation.
        """
        self.noise = []
        self.samplePoint = samplePoint
        self.num = endTrace - startTrace
        for i in range(DES_PI.NK):
            self.noise.append([])

        if (excludeStartTrace ==None and excludeEndTrace == None):
            exclude = False
        else:
            exclude = True

        if (endTrace == 0):
            endTrace = self.acq.numTraces

        for i in range(endTrace - startTrace):
            if (exclude == True and i>= excludeStartTrace and i<=excludeEndTrace):
                pass
            else:
                ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
                ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
                sboxIn = ptBlock ^ subkey
                sboxOut = des_block.sbox(sboxNum, sboxIn)
                hyp = HW(sboxOut)
                self.noise[hyp].append(self.acq.traces[i + startTrace][samplePoint])


        for i in range(DES_PI.NK):
            self.mean[i] = np.mean(self.noise[i])
            self.sd[i] = np.std(self.noise[i])
            #print "Num Measurements for sbox input %d is %d" % (i, len(noise[i]))
            if len(self.noise[i]) == 0:
                print "Not enough measurements for sbox input %d ! Please rebuild the template with more measurements" % i

        print self.mean
        print self.sd

    def ptNumpyArray2String(self, inputtext):
        ptStr = ""
        for i in range(self.acq.blockSize):
            ptStr += hex2str(inputtext[i])
        return ptStr

    def plot(self):
        yV = []

        xV=[]
        for i in range(0, self.num):
            xV.append(self.acq.traces[i][self.samplePoint])
        xV = set(xV)
        xV=list(xV)
        xV.sort()

        for i in range(5):
            yV.append([])
        #xV = range( min(x1,x2), max(x1,x2))
        for x in xV:
            for i in range(5):
                yV[i].append(scipy.stats.norm(self.mean[i], self.sd[i]).pdf(x))
        for i in range(5):
            plt.plot(xV, yV[i])
        plt.show()



if __name__ == "__main__":
    correctKey = "DEADBEEFCAFEBABE"
    sboxNum = 1
    roundNum = 1    # Do not Change this. Currently fixed to 1.
    rk = des_block.roundkey(correctKey,roundNum)   # Correct Round Key of roundNum
    sk = des_block.subkey(rk,sboxNum)             # Correct Sbox Key of roundNum, sboxNum
    samplePoint = 0

    acq = LoadTraces.LoadTraces('traces/0001.trs', samplePoint = 3881, numTraces = 100000)
    t = DES_PI(acq)
    #acq.shuffle()
    t.buildTemplate(0x1f, sboxNum, samplePoint, startTrace = 0, endTrace =100000)
    print "PI is %f" % (t.pi(0x1f, sboxNum, samplePoint, startTrace = 0, endTrace = 100000))


