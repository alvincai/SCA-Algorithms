#! /usr/bin/python

import LoadTraces
import des_block
import numpy as np
import matplotlib.pyplot as plt
from util import *
import scipy.stats
from collections import Counter
import pickle



def mod(old):
    newL = []
    for i in old:
        new = i - 2.3219280948873622
        new/= 16
        new *= 5
        new += 2.3219280948873622
        newL.append(new)
    return newL

class DES_PI():
    """ Calculate the percieved information for DES """

    NK = 5     # Number of DES Keys for an sbox
    H_x = 2.3219280948873622     # Entropy of sbox output
    Pr_x = 1.0/ 5.0    # Probability of an sbox output

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

        for i in range(startTrace, endTrace):
            norm = 0.0
            for sboxCandidate in range(DES_PI.NK):
                # normalise the sum of probabilities to 1
                norm += scipy.stats.norm(self.mean[sboxCandidate], self.sd[sboxCandidate]).pdf(self.acq.traces[i][samplePoint])

            for sboxIn in range(DES_PI.NK):
                PrXL = scipy.stats.norm(self.mean[sboxIn], self.sd[sboxIn]).pdf(self.acq.traces[i][samplePoint]) / norm  #Pr[X| Leakage] is estimated with this probabilty
                PrLX = PrXL * 5.0 / Nm

                if(PrXL!=0):
                    pi += DES_PI.Pr_x * PrLX * np.log2(PrXL)
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
                hyp = HD(sboxOut,des_block.ptRoundIn(ptStr,sboxNum))
                self.noise[hyp].append(self.acq.traces[i + startTrace][samplePoint])


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
    trsFile = 'traces/1001.trs'
    correctKey = "DEADBEEFCAFEBABE"
    sboxNum = 1

    # init
    roundNum = 1    # Do not Change this. Currently fixed to 1.
    rk = des_block.roundkey(correctKey,roundNum)   # Correct Round Key of roundNum
    sk = des_block.subkey(rk,sboxNum)             # Correct Sbox Key of roundNum, sboxNum
    acq = LoadTraces.LoadTraces(trsFile, numTraces=340001, samplePoint = 282)

    # Mutual Information
    pi = DES_PI(acq)     # init class
    xPI = []
    yPI = []
    xPImean = []
    yPImean = []

    numTraces = 50000
    piInterval = 1000
    for i in range(1000, numTraces, piInterval):
        meany = 0
        count = 0.0
        smallInterval = int(np.floor((i / 10.0)))
        for j in range(1,10):
            start = (j * smallInterval)
            end = (smallInterval*(j+1))
            pi.buildTemplate(sk, sboxNum, 0, startTrace = 0, endTrace=i, excludeStartTrace=start, excludeEndTrace = end )
            xPI.append(i)
            yval = pi.pi(sk, sboxNum, 0, startTrace = start, endTrace = end)
            yPI.append(yval)
            meany += yval
            count += 1.0
        meany /= count
        yPImean.append(meany)
        xPImean.append(i)
        acq.shuffle()

    pickle.dump(yPI, open('results/1001/yPI-HDRound.p', 'wb'))
    pickle.dump(xPI, open('results/1001/xPI-HDRound.p', 'wb'))
    pickle.dump(yPImean, open('results/1001/yPImean-HDRound.p', 'wb'))
    pickle.dump(xPImean, open('results/1001/xPImean-HDRound.p', 'wb'))

    #figure('Perceived Information')
    #plot(xPI, yPI, 'ro')
    #plot(xPImean, yPImean)
    #ylabel('Perceivec Information')
    #xlabel('Number of Measurements')
    #axisNew = []
    #axisOld = axis()
    #axisNew.append(0.0)
    #axisNew.append(numTraces + numTraces/10)
    #axisNew.append(axisOld[2])
    #axisNew.append(axisOld[3])
    #axis(axisNew)
    #draw()

    print "PI for 1001 is %f " % meany

    trsFile = 'traces/1002.trs'
    correctKey = "DEADBEEFCAFEBABE"
    sboxNum = 1
    leakage = 'HDRound'

    # init
    trsFile = 'traces/1002.trs'
    roundNum = 1    # Do not Change this. Currently fixed to 1.
    rk = des_block.roundkey(correctKey,roundNum)   # Correct Round Key of roundNum
    sk = des_block.subkey(rk,sboxNum)             # Correct Sbox Key of roundNum, sboxNum
    acq = LoadTraces.LoadTraces(trsFile, samplePoint=44, numTraces=340001)

    # Mutual Information
    pi = DES_PI(acq)     # init class
    xPI = []
    yPI = []
    xPImean = []
    yPImean = []

    numTraces = 50000
    piInterval = 1000
    for i in range(1000, numTraces, piInterval):
        meany = 0
        count = 0.0
        smallInterval = int(np.floor((i / 10.0)))
        for j in range(1,10):
            start = (j * smallInterval)
            end = (smallInterval*(j+1))
            pi.buildTemplate(sk, sboxNum, 0, startTrace = 0, endTrace=i, excludeStartTrace=start, excludeEndTrace = end )
            xPI.append(i)
            yval = pi.pi(sk, sboxNum, 0, startTrace = start, endTrace = end)
            yPI.append(yval)
            meany += yval
            count += 1.0
        meany /= count
        yPImean.append(meany)
        xPImean.append(i)
        acq.shuffle()
    pickle.dump(yPI, open('results/1002/yPI-HDRound.p', 'wb'))
    pickle.dump(xPI, open('results/1002/xPI-HDRound.p', 'wb'))
    pickle.dump(yPImean, open('results/1002/yPImean-HDRound.p', 'wb'))
    pickle.dump(xPImean, open('results/1002/xPImean-HDRound.p', 'wb'))

    print "PI for 1002 is %f " % meany
