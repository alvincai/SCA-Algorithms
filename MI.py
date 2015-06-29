#! /usr/bin/python

import LoadTraces
import des_block
import numpy as np
import matplotlib.pyplot as plt
from util import *
import scipy.stats
from collections import Counter

class DES_PI():
    """ Calculate the percieved information for DES """

    NK = 64     # Number of DES Keys for an sbox
    H_x = 6.0     # Entropy of sbox input (x)
    Pr_x = 1.0/ 64.0    # Probability of an sbox input

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

        leakage = []
        for i in range(startTrace, endTrace):
            leakage.append(self.acq.traces[i][samplePoint])
        leakage = set(leakage)

        for i in leakage:
            norm = 0.0
            for sboxCandidate in range(DES_PI.NK):
                # norm is the sum of Pr[L|X] for all X
                norm += scipy.stats.norm(self.mean[sboxCandidate], self.sd[sboxCandidate]).pdf(i)

            for sboxIn in range(DES_PI.NK):
                PrLX = scipy.stats.norm(self.mean[sboxIn], self.sd[sboxIn]).pdf(i) # Pr[L|X]
                PrXL = PrLX / norm  #Pr[X|L]

                if(PrXL!=0):
                    pi += DES_PI.Pr_x * PrLX * np.log2(PrXL)
        return pi

    def check(self):
        for i in self.mean:
            if np.isnan(i):
                print "Error, not enough values collected"

    def maxlikelihood(self, sboxNum, samplePoint, startTrace=0, endTrace=0):
        """
        Use the generated templates in a template attack to determine most likely key.
        Not related to PI (just for fun)
        First call buildTemplate to generate means and covariance matrix that will be used in this attack

        # Inputs:
        sboxNum - the sbox number we are targeting (values 1 to 8)
        samplePoint - the single sample point to calculate the PI for
        startTrace, endTrace - SNR is calculated from Trace startTrace to endTrace. This value is optional. If left out, it will calculate SNR for entire trace set

        Returns distinguisher values of sbox keys using Template attack, maximum likelihood
        """
        prob = np.ones(64, dtype=np.float64)

        for i in range(startTrace, endTrace):
            for sboxIn in range(DES_PI.NK):
                noise = self.acq.traces[i][samplePoint] - self.mean[sboxIn]
                p = scipy.stats.norm.pdf(noise, 0, self.sd[sboxIn]) # input to pdf are x, mean, sigma

                # from sboxinput and plaintext, decide which subkey the possibility belongs to
                ptStr = self.ptNumpyArray2String(self.acq.inputtext[i])
                ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
                sk = sboxIn ^ ptBlock
                prob[sk] *= p

        print "key is %x" % (np.argmax(prob))
        return prob



    def buildTemplate(self, subkey, sboxNum, samplePoint, startTrace = 0, endTrace = 0, excludeStartTrace =None, excludeEndTrace=None):
        """
        Calculate template for each sbox input
        a. Calculate mean of the signal
        b. Calculate variance of the noise

        # Inputs:
        subkey - the correct subkey
        sboxNum - the sbox number we are targeting (values 1 to 8)
        samplePoint - the single sample point to calculate the PI for
        startTrace, endTrace - SNR is calculated from Trace startTrace to endTrace. This value is optional. If left out, it will calculate SNR for entire trace set
        excludestartTrace, excludedEndTrace - Exclude the traces between excludeStartTrace to excludeEndTrace. This interval should be between startTrace and endTrace. It is used for cross validation.
        """
        self.noise = []
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
                self.noise[sboxIn].append(self.acq.traces[i + startTrace][samplePoint])



        for i in range(DES_PI.NK):
            self.mean[i] = np.mean(self.noise[i])
            self.sd[i] = np.std(self.noise[i])
            #print "Num Measurements for sbox input %d is %d" % (i, len(noise[i]))
            if len(self.noise[i]) == 0:
                print "Not enough measurements for sbox input %d ! Please rebuild the template with more measurements" % i



    def ptNumpyArray2String(self, inputtext):
        ptStr = ""
        for i in range(self.acq.blockSize):
            ptStr += hex2str(inputtext[i])
        return ptStr



if __name__ == "__main__":
    acq = LoadTraces.LoadTraces('../thesis/traces/0002.trs', samplePoint = 18694)
    t = DES_PI(acq)
    sboxNum = 1
    samplePoint = 0

    t.buildTemplate(0x1f, sboxNum, samplePoint, startTrace = 0, endTrace =9000)
    print "PI is %f" % (t.pi(0x1f, sboxNum, samplePoint, startTrace = 9501, endTrace = 10000))


