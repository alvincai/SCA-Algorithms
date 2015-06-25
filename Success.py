#! /usr/bin/python

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
import timer

class DESSuccess():

    """Calculate the Rank0 success rate, both empirical and theoretical (using Confusion analysis)"""

    # Lower limit of the integration for calculating cdf (used in scipy.stats.mvn)
    lower = np.zeros(63)
    # Upper limit of the integration for calculating cdf (used in scipy.stats.mvn)
    upper = np.ones(63) * 5000



    def __init__(self, acq, sboxNum, correctSubKey, signal, noise, leakage='HW'):
        self.acq = acq
        self.sboxNum = sboxNum
        self.correctSubKey = correctSubKey
        self.signal = signal
        self.noise = noise
        self.leakage = leakage
        self.dpa = DES_DPA.DES_DPA(acq)
        self.cpa = DES_CPA.DES_CPA(acq)
        self.conf = Confusion.DESConfusion(sboxNum, leakage = leakage)



    def TheoreticalDPA(self, numTraces):
        """
        Returns the DES, DPA Success Rate (calculated using Yunsi Fei's Confusion Analysis formula) for the indicated number of traces
        """
        cov, mean = self.conf.YCovMeanDPA(self.signal, self.noise, numTraces, correctKey=self.correctSubKey)
        return mvn.mvnun(self.lower, self.upper, mean.reshape(63), cov)[0]


    def TheoreticalCPA(self, numTraces):
        """
        Returns the DES, CPA Success Rate (calculated using Yunsi Fei's Confusion Analysis formula) for the indicated number of traces
        """
        cov, mean = self.conf.YCovMeanCPA2(self.signal, self.noise, numTraces, correctKey=self.correctSubKey)
        return mvn.mvnun(self.lower, self.upper, mean.reshape(63), cov)[0]


    def Empirical(self, numTraces, samplePoint = None, attack='cpa', runs=1000, leakage='HW'):
        """
        Returns Empirical success rate for DPA for the indicated Number of measurements (numTraces)
        If samplePoint is provided, only uses single samplePoint, else use entire trace to calculate distinguishers
        Performs attack 1000 times and returns the number of rank0 success/1000

        # inputs:
        numTraces - number of traces to use to calculate empirical success rate
        samplePoint -  the single samplePoint of the leakage which we use to calculate the success rate
        attack - input 'cpa' or 'dpa' depending on the desired distinguisher type
        runs - how many times to conduct the attack. Typically use 1000 or 100 for better statistical properties
        leakage - The leakage model. This takes values 'HW', 'HD', 'HDRound'

        # Output:
        Empirical sucess rate as probability (between 0 and 1)
        """
        success = 0
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
            maxAbs = np.argmax(np.absolute(distg))
            if (self.correctSubKey == maxAbs):
                # success
                success += 1.0

        return success/runs

    def verify(self):
        """ Verify our results against Yunsi Fei's results where she publishes her signal to noise ration and success rate curves
        """
        x = []
        yDPA = []
        yCPA = []

        # Amended Formula which works
        for i in range(1,200,20):
            cov, mean = self.conf.YCovMeanDPA(0.0016, 0.0046, i, correctKey=60)
            yDPA.append(mvn.mvnun(DESSuccess.lower, DESSuccess.upper, mean.reshape(63), cov)[0])
            x.append(i)
            cov, mean = self.conf.YCovMeanCPA2(0.0016, 0.0048, i, correctKey=60)
            yCPA.append(mvn.mvnun(self.lower, self.upper, mean.reshape(63), cov)[0])

        p1 = plot(xest,ySR1est,'--')
        p2 = plot(x,yCPA)
        legend(loc='lower right')
        ylabel('Success Rate')
        xlabel('Number of Measurements')
        show()

        # Paper's formula which does not work!!
        #lower2 = np.ones(63) * (-1000)
        #meanNorm = np.zeros(63)

        #for i in range(10,1000,10):
            #cov, mean = self.conf.YCovMean(self.signal, self.noise, i, correctKey=60)
            ## Need to test if the inv(sqrtm(cov)) computation is computed correctly.
            ## If it is not, we discard the result as it is meaningless
            #if (self.test(cov)):
                #val = (i ** 0.5) * np.dot(mean, inv(sqrtm(cov)))
                #y.append(mvn.mvnun(lower2, val,  meanNorm, np.eye(63)))
                #x.append(i)

        #plt.plot(x,y)
        #plt.ylabel('Success Rate')
        #plt.xlabel('Number of Measurements')
        #plt.show()


    def test(self, cov):
        """
        Test if the inv(sqrt(cov)) can be calculated correctly by numpy
        This uses the property that A x inv(A) = identity matrix
        Returns the result of this test
        """
        A = sqrtm(cov)
        invA = inv(A)
        return (np.allclose(np.eye(63), np.dot(A, invA), rtol = 0.1, atol = 0.1))

    def Cdf(self,x):
        # cdf = integral of pdf from lower limit to upper limit (x)
        # Reference http://www.nhsilbert.net/source/2014/04/multivariate-normal-cdf-values-in-python/
        return mvn.mvnun(self.lowerLimit63, x, self.meanNorm63, self.covNorm63)

    def lsb(self, num):
        return num & 1

    def msb(self, num):
        return num & 0x8


if __name__ == "__main__":

    # Example for calculating the theoretical success rate
    acq = None  #if calculating theoretical and not emprical success rate, just feed None into acq
    sboxNum = 1
    correctSubKey = 0x00
    signal = 4.0
    noise = 22.85
    numTraces = 1500 # the number of traces for the theoretical success rate plot
    interval = numTraces/20
    s1 = DESSuccess(acq, sboxNum, correctSubKey, signal, noise, leakage='HW')


    ySR = []
    xSR = []
    for i in range (1, numTraces, interval):
        ySR.append(s1.TheoreticalCPA(i))
        xSR.append(i)

    figure('Success Rate')
    plot(xSR, ySR)
    xlabel('Number of Measurements')
    ylabel('')
    show()



