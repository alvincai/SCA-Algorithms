#! /usr/bin/python
# Based on the papers:
# [1] Fei, Yunsi, Qiasi Luo and Adam Ding. "A statistical model for DPA with novel algorithmic confusion analysis." Cryptographic Hardware and Embedded Systems-CHES2012. Springer Berlin Heidelberg, 2012. 233-250."
# [2] Fei, Yunsi, et al. "A Statistics-based Fundamental Model for Side-channel Attack Analysis." IACR Cryptology ePrint Archive 2014 (2014): 152.


import des_block
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
from util import *

class DESConfusion():
    """
    This Class generates DES Confusion coefficients necessary for success rate simulation.
    It is not meant to be used in a standalone manner.
    Instead the results from this class are used by the class, Success.
    """


    def __init__(self, sboxNum, leakage='HW'):
        self.NkDES = 2**6
        self.sboxNum = sboxNum
        self.leakage = leakage
        self.CPAConfusionCoefficients = self.DESCpaCC(sboxNum)  # CPA Confusion Co-efficients
        self.DPAConfusionCoefficients = self.DESDpaCC(sboxNum)  # DPA Confusion Co-efficients
        self.dpaCV = self.DPA_CV()  # DPA Confusion Vector for key = 0
        self.dpaK = self.DPA_K()    # DPA Confusion Matrix for key = 0
        self.cpaCV = self.CPA_CV()  # CPA Confusion Vector
        self.cpaK, self.cpaKs, self.cpaKss = self.CPA_K()    # CPA Confusion Matrix, K, K* and K**

    def YCovMeanDPA(self, signal, noise, Nm, correctKey=0):
        """
        This function calculates the mean and covariance matrix of the multivariate gaussian variable Y
        This is used directly in the Success class to calculate SR for DES DPA
        Refer to 'Statistical model for DPA with Confusion Analysis'
        # input:
          Signal and noise values
          Number of measurements (Nm)
          correctKey: (0...63) -- Note: Correct key does not matter because of symmetry
        # Return :
          covariance matrix of size [(Nk-1) x (Nk-1)]
          mean of size [(Nk-1) x 1]
        """

        Nk = self.NkDES
        if (correctKey != 0):
            cv = self.DPA_CV(correctKey = correctKey)  # Confusion Vector
            K = self.DPA_K(correctKey = correctKey)    # Confusion Matrix
        else:
            cv = self.dpaCV
            K = self.dpaK

        cov = np.zeros((Nk-1, Nk-1))
        mean = cv * 2 * signal
        for i in range(Nk-1):
            for j in range(Nk-1):
                ki = cv[0][i]
                kj = cv[0][j]
                cov[i][j] = ((16 * K[i][j] * (noise**2)) + (4 * (K[i][j] - (ki * kj)) *(signal**2)))/ Nm

        return cov, mean

    def YCovMeanCPA2(self, signal, noise, Nm, correctKey=0):
        """
        This method calculates the mean and covariance matrix of the multivariate guassian variable Y
        This is used directly in the Success class to calculate SR for DES CPA
        Refer to 'A Statistics-based Fundamental Model for Side-channel Attack Analysis'
        # input:
          Signal and noise values
          Number of measurements (Nm)
          correctKey: (0...63) -- Note: Correct key does not matter because of symmetry
        # Return :
          covariance matrix of size [(Nk-1) x (Nk-1)]
          mean of size [(Nk-1) x 1]
        """

        Nk = self.NkDES
        if (correctKey != 0):
            cv = self.CPA_CV(correctKey = correctKey)       # Confusion Vector
            K, Ks, Kss = self.CPA_K(correctKey = correctKey)     # Confusion Matrices K, K**
        else:
            cv = self.cpaCV
            K = self.cpaK
            Ks = self.cpaKs
            Kss = self.cpaKss
        cov = np.zeros((Nk-1, Nk-1))
        snr = signal / noise
        mean = 0.5 * (snr **2 ) * cv

        for i in range(Nk-1):
            for j in range(Nk-1):
                ki = cv[0][i]
                kj = cv[0][j]
                #cov[i][j] = ((K[i][j] * (snr**2)) )/ Nm    # Simplified formula for low SNR
                cov[i][j] = ((K[i][j] * (snr**2)) + (0.25 * (Kss[i][j] - (ki * kj)) *(snr ** 4)))/ Nm
        return cov, mean

    def DESDpaCC(self, sboxNum, plot = False):
        """ Returns DES DPA Confusion coefficents Matrix of size [number of keys x number of keys]
        #inputs:
        sboxNum - the sbox number we are targeting (values 1 to 8)
        plot - if True, will plot the matrix as a histogram
        """
        numPT = 2 ** 6      # number of possible PT/CT
        numKey = 2 ** 6     # number of possible keys
        cc = np.zeros((numKey,numKey), np.float)    # Confusion Coefficient matrix
        histogram = []

        for ki in range(numKey):
            for kj in range(numKey):
                numNotEqual = 0
                for ptBlock in range(numPT):
                    sboxOuti = des_block.sbox(sboxNum, ptBlock ^ ki)
                    sboxOutj = des_block.sbox(sboxNum, ptBlock ^ kj)
                    if (self.msb(sboxOuti) != self.msb(sboxOutj)):
                        numNotEqual += 1.0
                coefficient = numNotEqual / numPT
                cc[ki][kj] = coefficient
                if (ki != kj and ki<kj):
                    histogram.append(coefficient)

        if (plot):
            # Plot a histogram of the coefficients
            weights = np.ones_like(histogram)/len(histogram)
            fig = plt.hist(histogram, 1000, weights=weights)
            plt.ylabel('Frequency')
            plt.xlabel('Confusion coefficient')
            plt.show(fig)
        return cc

    def DESCpaCC(self, sboxNum, plot = False, leakage ='HW'):
        """ Returns DES CPA 2-way Confusion coefficents Matrix of size [number of keys x number of keys]
        #inputs:
        sboxNum - the sbox number we are targeting (values 1 to 8)
        plot - if True, will plot the matrix as a histogram
        leakage - The leakage model. This takes values 'HW', 'HD', 'HDRound'
        """

        numPT = 2 ** 6      # number of possible PT/CT
        numKey = 2 ** 6     # number of possible keys
        cc = np.zeros((numKey,numKey), np.float)    # Confusion Coefficient matrix
        histogram = []

        for ki in range(numKey):
            for kj in range(numKey):
                numNotEqual = 0.0
                k = []
                for ptBlock in range(numPT):
                    sboxIni = ptBlock ^ ki
                    sboxInj = ptBlock ^ kj

                    sboxOuti = des_block.sbox(sboxNum, sboxIni)
                    sboxOutj = des_block.sbox(sboxNum, sboxInj)
                    if leakage =='HW':
                        k.append((self.hw(sboxOuti) - self.hw(sboxOutj)) ** 2)
                    if leakage =='HD':
                        k.append((HD(sboxOuti, sboxIni) - HD(sboxOutj,sboxInj)) ** 2)


                cc[ki][kj] = np.mean(k)
                if (ki != kj and ki<kj):
                    histogram.append(cc[ki][kj])

        if (plot):
            weights = np.ones_like(histogram)/len(histogram)
            fig = plt.hist(histogram, 1000, weights=weights)
            plt.ylabel('Frequency')
            plt.xlabel('Confusion coefficient')
            plt.show(fig)

        return cc

    def DPA_CV(self, correctKey=0):
        """ Return the DES DPA confusion vector of the correct Key. This is a Nk x1 matrix """
        Nk = self.NkDES
        cv = np.delete(self.DPAConfusionCoefficients[correctKey], correctKey)
        cv = cv.reshape(1, Nk-1)
        return cv

    def CPA_CV(self, correctKey=0):
        """ Return the DES CPA confusion vector of the correct Key. This is a Nk x1 matrix """
        Nk = self.NkDES
        cv = np.delete(self.CPAConfusionCoefficients[correctKey], correctKey)
        cv = cv.reshape(1, Nk-1)
        return cv

    def DPA_K(self, correctKey=0):
        """ Returns the Nk-1 x Nk-1 Confusion matrix
        #input:
        correctKey - takes values (0...63)
        """
        Nk = self.NkDES
        cc2 = np.delete(self.DPAConfusionCoefficients, correctKey,0)
        cc2 = np.delete(cc2, correctKey,1)      # Tmp variable used to calculate the 3-way co-efficient
        cv = self.DPA_CV(correctKey = correctKey)
        K = np.zeros((Nk-1, Nk-1))

        for i in range(Nk-1):
            for j in range(Nk-1):
                if (i == j):
                    K[i][j] = cv[0][i]
                else:
                    ki = cv[0][i]
                    kj = cv[0][j]
                    kij = cc2[i][j]
                    kcij = 0.5 * (ki + kj - kij)    # 3 way confusion co-efficient
                    K[i][j] = kcij
        return K

    def CPA_K(self, correctKey=0, leakage='HW'):
        """Returns the 3-way CPA confusion matrix K, K* and K**
        (refer to "A Statistics-based Fundamental Model for Side-channel Attack Analysis.")

        #inputs:
        correctKey - takes values (0...63)
        leakage - the leakage model. Takes values 'HW' and 'HD'.
        """
        Nk = self.NkDES
        numPT = 2 ** 6
        sboxNum = self.sboxNum
        cv = self.CPA_CV(correctKey = correctKey) # Diagonal of the confusion matrix
        K = np.zeros((Nk, Nk))
        Ks = np.zeros((Nk, Nk))     # K*
        Kss = np.zeros((Nk, Nk))    # K**
        keys = np.arange(64)    # List of wrong keys
        keys = np.delete(keys, correctKey)

        evkc = []    # E [V | kc]
        for ptBlock in range(numPT):
            sboxOutc = des_block.sbox(sboxNum, ptBlock ^ correctKey)
            evkc.append(self.hw(sboxOutc))
        evkc = np.mean(evkc)
        #print evkc

        for i in keys:
            for j in keys:
                # Calculate kcij = E[(V|kc - V|ki) * (V|kc - V|ki)]
                # Calculate kcijss = E[4 * (V|kc - E[V|kc])^2 * (V|kc - V|ki) * (V|kc - V | kj)]
                kcij = []
                kcijs = []
                kcijss = []
                for ptBlock in range(numPT):
                    sboxOutc = des_block.sbox(sboxNum, ptBlock ^ correctKey)
                    sboxOuti = des_block.sbox(sboxNum, ptBlock ^ i)
                    sboxOutj = des_block.sbox(sboxNum, ptBlock ^ j)
                    if self.leakage =='HW':
                        vkc = self.hw(sboxOutc) # V | kc
                        vki = self.hw(sboxOuti) # V | ki
                        vkj = self.hw(sboxOutj) # V | kj
                    elif self.leakage =='HD':
                        vkc = HD(sboxOutc, correctKey^ptBlock) # V | kc
                        vki = HD(sboxOuti, i^ptBlock) # V | ki
                        vkj = HD(sboxOutj, j^ptBlock) # V | kj


                    kcij.append((vkc - vki) * (vkc -vkj))
                    kcijs.append(((vkc-vki)**2) * ((vkc-vkj)**2))
                    kcijss.append( 4 * ((vkc - evkc)**2) * (vkc - vki) * (vkc -vkj))
                kcij = np.mean(kcij)
                kcijss = np.mean(kcijss)
                kcijs = np.mean(kcijs)

                K[i][j] = kcij
                Ks[i][j] = kcijs
                Kss[i][j] = kcijs

        K = np.delete(K, correctKey,0)
        K = np.delete(K, correctKey,1)
        Ks = np.delete(Ks, correctKey,0)
        Ks = np.delete(Ks, correctKey,1)
        Kss = np.delete(Kss, correctKey,0)
        Kss = np.delete(Kss, correctKey,1)

        return K, Ks, Kss

    def hw(self, num):
        return bin(num).count("1")

    def lsb(self, num):
        return num & 1

    def msb(self, num):
        return num & 0x8

    def DesCPA(self, sboxNum):
        pass

    def is_pos_def(self,covx):
        return np.all(np.linalg.eigvals(cov) > 0)

if __name__ == "__main__":
    sboxNum = 1
    c = DESConfusion(sboxNum)
    c.DESCpaCC(1, plot= True)
    #cov, mean = c.YCovMeanDPA(0.0016, 0.0048, 10, correctKey=0)


