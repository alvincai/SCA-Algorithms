#! /usr/bin/python
import LoadTraces
import des_block
import DES_CPA
import DES_DPA
import numpy as np
import matplotlib.pyplot as plt
from util import *
from matplotlib.pyplot import *
import pickle


class SNR():
    """This class calculates the CPA or DPA Signal to noise ratio
    Note: Only round 1 of DES
    """

    def __init__(self, acq):
        self.acq = acq

    def cpaTrace(self, subkey, sboxNum, startTrace =0, endTrace = 0, leakage='HW'):
        """The cpaTrace function calculates the CPA signal to noise ratio over all the samples of a traceset
        Inputs:
        subkey - the correct subkey (which is 0 to 63 for DES)
        sboxNum - the sbox number we are targeting (values 1 to 8)
        startTrace, endTrace - SNR is calculated from Trace startTrace to endTrace. This value is optional. If left out, it will calculate SNR for entire trace set
        leakage - The leakage model. This takes values 'HW', 'HD', 'HDRound'

        Returns:
        Signal matrix, noise matrix
        """

        sumnum = np.zeros(self.acq.numSamples)
        if (endTrace == 0):
            endTrace = self.acq.numTraces


        # Group traces according to the hamming weight of the leakage
        traceHW = []
        noise = []
        for i in range(5):
            traceHW.append(None)
            noise.append(None)

        # Generate hypotheticals
        hyp = np.zeros(endTrace - startTrace, dtype=int)
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

            if traceHW[hyp[i]] == None:
                traceHW[hyp[i]] = self.acq.traces[i]
            else:
                traceHW[hyp[i]] = np.vstack((traceHW[hyp[i]] ,self.acq.traces[i]))


        meanh = np.mean(hyp, dtype = np.float)
        meant = np.mean(self.acq.traces[startTrace:endTrace], axis=0, dtype=np.float)
        noise = np.std(traceHW[0], axis=0)
        for i in range(5):
            noise = np.vstack((noise, np.std(traceHW[i], axis=0)))
        noise = np.mean(noise, axis = 0)


        for tnum in range(endTrace - startTrace):
            hdiff = (hyp[tnum] - meanh)
            tdiff = self.acq.traces[tnum+startTrace,:] - meant

            sumnum += hdiff * tdiff
            signal = np.abs(sumnum / (endTrace - startTrace))

        return signal, noise

    def cpaPoint(self, sk, sboxNum, samplePoint, startTrace = 0, endTrace=0, leakage='HW'):
        """The cpaPoint function calculates the CPA signal to noise ratio for a single sample of a given traceset.
        This saves a lot of time compared to cpaTrace.
        The recommended usage scenario will be to determine the sample point of interest using cpaTrace (with fewer traces),
        and then to use cpaPoint (with many traces)

        Inputs:
        subkey - the correct subkey (which is 0 to 63 for DES)
        sboxNum - the sbox number we are targeting (values 1 to 8)
        samplePoint - the single samplePoint to calculate the SNR for
        startTrace, endTrace - SNR is calculated from Trace startTrace to endTrace. This value is optional. If left out, it will calculate SNR for entire trace set
        leakage - The leakage model. This takes values 'HW', 'HD', 'HDRound'

        returns:
        signal, noise
        """
        if endTrace == 0:
            endTrace = self.acq.numTraces

        HammW = []
        for i in range(5):
            HammW.append([])

        numTraces = endTrace - startTrace
        covMethod = np.zeros((2,numTraces))

        for i in range(numTraces):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxOut = des_block.sbox(sboxNum, ptBlock ^ sk)

            if leakage =='HW':
                hyp = HW(sboxOut)   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HDRound':
                hyp = HD(sboxOut,des_block.ptRoundIn(ptStr,sboxNum))   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HD':
                hyp  = HD(sboxOut, sboxIn)
            else:
                print "undefined leakage model"
                return 0

            HammW[hyp].append(self.acq.traces[i+startTrace][samplePoint])
            covMethod[0][i] = self.acq.traces[i+startTrace][samplePoint]
            covMethod[1][i] = hyp

        noiseCPA = []
        signalCPA = []
        signalCPAdiff = []

        for j in range(5):
            noise = np.std(HammW[j])
            if (np.isnan(noise)):
                pass
            else:
                noiseCPA.append(noise)
        noiseCPA = np.mean(noiseCPA)
        signalCov = np.abs(((np.cov(covMethod))[1][0]))

        return signalCov*self.acq.yscale, noiseCPA *self.acq.yscale

    def cpaPoint2(self, sk, sboxNum, samplePoint, startTrace = 0, endTrace=0, leakage='HW'):
        """Same as cpaPoint but uses the traditional `means' method instead of the covariance method.  """
        if endTrace == 0:
            endTrace = self.acq.numTraces

        HammW = []
        for i in range(5):
            HammW.append([])

        numTraces = endTrace - startTrace
        #covMethod = np.zeros((2,numTraces))

        for i in range(numTraces):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxOut = des_block.sbox(sboxNum, ptBlock ^ sk)

            if leakage =='HW':
                hyp = HW(sboxOut)   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HDRound':
                hyp = HD(sboxOut,des_block.ptRoundIn(ptStr,sboxNum))   # or use HD(sboxOut, sboxIn) for hardware implementations
            elif leakage =='HD':
                hyp  = HD(sboxOut, sboxIn)
            else:
                print "undefined leakage model"
                return 0

            HammW[hyp].append(self.acq.traces[i+startTrace][samplePoint])

        noiseCPA = []
        signalCPA = []
        signalCPAdiff = []

        for j in range(5):
            noise = np.std(HammW[j])
            if (np.isnan(noise)):
                print "noise isnan, rejected"
                pass
            else:
                noiseCPA.append(noise)
            signalCPA.append(np.mean(HammW[j]))
        noiseCPA = np.mean(noiseCPA)
        for j in range(4):
            signalCPAdiff.append(np.abs(signalCPA[j+1] - signalCPA[j]))
        signalCPA = np.mean(signalCPAdiff)
        return signalCPA*self.acq.yscale, noiseCPA*self.acq.yscale


    def dpaPoint(self, sk, sboxNum, samplePoint, startTrace = 0, endTrace=0, bit=1):
        """The dpaPoint function calculates the DPA signal to noise ratio for a single sample of a given traceset.

        Inputs:
        subkey - the correct subkey (which is 0 to 63 for DES)
        sboxNum - the sbox number we are targeting (values 1 to 8)
        samplePoint - the single samplePoint to calculate the SNR for
        startTrace, endTrace - SNR is calculated from Trace startTrace to endTrace. This value is optional. If left out, it will calculate SNR for entire trace set
        leakage - The leakage model. This takes values 'HW', 'HD', 'HDRound'
        if bit is provided, returns the SNR of that specific sbox output bit. This takes values of 0,1,2,3
        """

        if endTrace == 0:
            endTrace = self.acq.numTraces


        LSB = []

        for i in range(2):
            LSB.append([])
        numTraces = endTrace - startTrace

        for i in range(numTraces):
            ptStr = self.ptNumpyArray2String(self.acq.inputtext[i + startTrace])
            ptBlock = des_block.ptblock2(ptStr, sboxNum)     # 6-bit block of plaintext which is fed into sboxNum
            sboxOut = des_block.sbox(sboxNum, ptBlock ^ sk)

            LSB[(sboxOut & (2**bit)) >> bit].append(self.acq.traces[i+startTrace][samplePoint])

        signalDPA = np.abs(np.mean(LSB[0]) - np.mean(LSB[1]))
        noiseDPA = []

        # Traditional method for calculating SNR for DPA
        for j in range(2):
            noiseDPA.append(np.std(LSB[j]))
        noiseDPA = np.mean(noiseDPA)

        return signalDPA * self.acq.yscale, noiseDPA * self.acq.yscale



    def ptNumpyArray2String(self, inputtext):
        """Format the inputtext into format acceptable to des_block"""

        ptStr = ""
        for i in range(self.acq.blockSize):
            ptStr += hex2str(inputtext[i])
        return ptStr

if __name__ == "__main__":
    trsFile = '/Data/Documents/0007.trs'
    acq = LoadTraces.LoadTraces(trsFile, samplePoint = 39205)
    sboxNum = 1
    roundNum = 1    # Do not change this
    correctKey = "DEADBEEFCAFEBABE"
    leakage = 'HW'  # Use HW, HD or HDRound
    rk = des_block.roundkey(correctKey,roundNum)   # Correct Round Key of roundNum
    sk = des_block.subkey(rk,sboxNum)             # Correct Sbox Key of roundNum, sboxNum

    s = SNR(acq)
    #sig, noi = s.cpaTrace(sk, sboxNum, startTrace = 0, endTrace=10000, leakage = leakage)
    #samplePoint = np.argmax(np.abs(sig/noi))

    #print "Signal is %f, Noise is %f" % (sig[samplePoint], noi[samplePoint])
    #print "SNR is %f" % (sig[samplePoint]/noi[samplePoint])



    # Verify the accuracy using cross validation
    # Delete redundant samples and also randomise the order of the traces
    acq.shuffle()

    numTraceSNR = 1000 # acq.numTraces/10
    xSNR = []
    ySignal = []
    yNoise = []
    ySNR = []
    ySNRMean = []
    xSNRMean = []
    yNoiseMean = []
    ySignalMean = []

    for i in range(numTraceSNR/20, numTraceSNR, numTraceSNR/10):
        meanSig = 0.0
        meanNoi = 0.0
        meanSNR = 0.0
        count = 0.0
        for j in range(1,10):
            signalTmp, noiseTmp = s.cpaPoint(sk, sboxNum, 0, startTrace=i*j, endTrace = (i*(j+1)), leakage=leakage)
            signalTmp *= 1000
            noiseTmp *= 1000
            xSNR.append(i)
            ySignal.append(signalTmp)
            yNoise.append(noiseTmp)
            ySNR.append(signalTmp/noiseTmp)
            meanSNR += signalTmp/noiseTmp
            meanSig += signalTmp
            meanNoi += noiseTmp
            count += 1.0
        meanSig /= count
        meanNoi /= count
        meanSNR /= count

        ySignalMean.append(meanSig)
        yNoiseMean.append(meanNoi)
        ySNRMean.append(meanSNR)
        xSNRMean.append(i)
        acq.shuffle()

    figure('Signal to Noise Ratio')
    ax = subplot(121)
    plot(xSNR, ySignal, 'ro', label='Signal')
    plot(xSNR, yNoise, 'b+', label='Noise')
    plot(xSNRMean, ySignalMean, 'r')
    plot(xSNRMean, yNoiseMean, 'b')
    ax.set_title('Signal and Noise', fontsize=16)
    legend(loc='upper right')
    ylabel('Signal and noise (mV)', fontsize=14)
    xlabel('Number of Measurements', fontsize=14)

    ax = subplot(122)
    ax.set_title('SNR', fontsize=16)
    plot(xSNR, ySNR, 'ro')
    plot(xSNRMean, ySNRMean, 'r')
    xlabel('Number of Measurements')
    show()


