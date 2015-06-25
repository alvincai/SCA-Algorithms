#! /usr/bin/python

import numpy as np
import Success
import SNR
import DES_DPA
import DES_CPA
import LoadTraces
import des_block
import RDM
import PI
from matplotlib.pyplot import *
import sys
import pickle
from timer import Timer

if __name__ == "__main__":

    # Fill in the appropriate tracefile name, Correct DES Key, DES round Number  and sbox number to target
    # Fill in the number of Traces. This is the upper limit for which the various metrics will be calculated. When calculating success rate, the number of traces in the trs file should be > numTraces * 1000
    # e.g. if numTraces = 100, for success rate, we will calculate success rate with numTraces =0 ... until success rate with numTraces = 100
    trsFile = 'traces/0001.trs'
    correctKey = "DEADBEEFCAFEBABE"
    sboxNum = 1
    numTraces = 30

    # init
    roundNum = 1    # Do not Change this. Currently fixed to 1.
    rk = des_block.roundkey(correctKey,roundNum)   # Correct Round Key of roundNum
    sk = des_block.subkey(rk,sboxNum)             # Correct Sbox Key of roundNum, sboxNum
    acq = LoadTraces.LoadTraces(trsFile, samplePoint = 3881, numTraces=1)
    #acq.shuffle()



    # Empirical Success Rate
    #runs = 20     # numberof attacks to do for empricial success rate calculation


    #sr = Success.DESSuccess(acq, sboxNum, sk, 1,1)
    #ySREmpirical = []
    #xSREmpirical = []

    #for i in range(5,4999,500):
        #with Timer() as t:
            #ySREmpirical.append(sr.Empirical(i, samplePoint=0, attack='cpa', runs=runs))
        #print "=> elasped sr for %d num of traces: %s s" % (i,t.secs)
        #ySREmpirical.append(t.secs*5)
        #xSREmpirical.append(i)

    #for i in range(5000,9999,1000):
        #with Timer() as t:
            #ySREmpirical.append(sr.Empirical(i, samplePoint=0, attack='cpa', runs=runs))
        #print "=> elasped sr for %d num of traces: %s s" % (i,t.secs)
        #ySREmpirical.append(t.secs*5)
        #xSREmpirical.append(i)

    #for i in range(10000,50001,10000):
        #with Timer() as t:
            #ySREmpirical.append(sr.Empirical(i, samplePoint=0, attack='cpa', runs=runs))
        #print "=> elasped sr for %d num of traces: %s s" % (i,t.secs)
        #ySREmpirical.append(t.secs*5)
        #xSREmpirical.append(i)

    #pickle.dump(ySREmpirical, open('results/ySREmpiricalTime.p', 'wb'))
    #pickle.dump(xSREmpirical, open('results/xSREmpiricalTime.p', 'wb'))




    ### Mutual Information
    #pi = PI.DES_PI(acq)     # init class
    #xPI = []
    #yPI = []

    #for i in range(1000, 4999, 500):
        #with Timer() as t:
            #smallInterval = int(np.floor((i / 10.0)))
            #j = 1
            #start = (j * smallInterval)
            #end = (smallInterval*(j+1))
            #pi.buildTemplate(sk, sboxNum, 0, startTrace = 0, endTrace=i, excludeStartTrace=start, excludeEndTrace = end )
            #yval = pi.pi(sk, sboxNum, 0, startTrace = start, endTrace = end)
        #print "=> elasped PI for %d num of traces: %s s" % (i,t.secs)
        #yPI.append(t.secs)
        #xPI.append(i)

    #for i in range(5000, 9999, 1000):
        #with Timer() as t:
            #smallInterval = int(np.floor((i / 10.0)))
            #j = 1
            #start = (j * smallInterval)
            #end = (smallInterval*(j+1))
            #pi.buildTemplate(sk, sboxNum, 0, startTrace = 0, endTrace=i, excludeStartTrace=start, excludeEndTrace = end )
            #yval = pi.pi(sk, sboxNum, 0, startTrace = start, endTrace = end)
        #print "=> elasped PI for %d num of traces: %s s" % (i,t.secs)
        #yPI.append(t.secs)
        #xPI.append(i)

    #for i in range(10000, 50001, 10000):
        #with Timer() as t:
            #smallInterval = int(np.floor((i / 10.0)))
            #j = 1
            #start = (j * smallInterval)
            #end = (smallInterval*(j+1))
            #pi.buildTemplate(sk, sboxNum, 0, startTrace = 0, endTrace=i, excludeStartTrace=start, excludeEndTrace = end )
            #yval = pi.pi(sk, sboxNum, 0, startTrace = start, endTrace = end)
        #print "=> elasped PI for %d num of traces: %s s" % (i,t.secs)
        #yPI.append(t.secs)
        #xPI.append(i)


    #pickle.dump(yPI, open('results/yPITime.p', 'wb'))
    #pickle.dump(xPI, open('results/xPITime.p', 'wb'))



    #SNR
    #s = SNR.SNR(acq)
    #ySNR =[]
    #xSNR=[]

    #for i in range(5,4999,500):
        #with Timer() as t:
            #s.cpaPoint(sk, sboxNum, 0, startTrace=0, endTrace = i, leakage='HW')
        #print "=> elasped SNR for %d num of traces: %s s" % (i,t.secs)
        #ySNR.append(t.secs)
        #xSNR.append(i)

    #for i in range(5000,9999,1000):
        #with Timer() as t:
            #s.cpaPoint(sk, sboxNum, 0, startTrace=0, endTrace = i, leakage='HW')
        #print "=> elasped SNR for %d num of traces: %s s" % (i,t.secs)
        #ySNR.append(t.secs)
        #xSNR.append(i)


    #for i in range(10000,50001,10000):
        #with Timer() as t:
            #s.cpaPoint(sk, sboxNum, 0, startTrace=0, endTrace = i, leakage='HW')
        #print "=> elasped SNR for %d num of traces: %s s" % (i,t.secs)
        #ySNR.append(t.secs)
        #xSNR.append(i)



    #pickle.dump(ySNR, open('results/ySNRTime.p', 'wb'))
    #pickle.dump(xSNR, open('results/xSNRTime.p', 'wb'))

    xSREmpirical = pickle.load(open('results/xSREmpiricalTime.p', 'rb'))
    ySREmpirical = pickle.load(open('results/ySREmpiricalTime2.p', 'rb'))
    xPI = pickle.load(open('results/xPITime.p', 'rb'))
    yPI = pickle.load(open('results/yPITime.p', 'rb'))
    xSNR = pickle.load(open('results/xSNRTime.p', 'rb'))
    ySNR = pickle.load(open('results/ySNRTime.p', 'rb'))

    matplotlib.rc('xtick', labelsize=15)
    matplotlib.rc('ytick', labelsize=15)





    figure('Time Taken')
    ax = subplot(131)
    p1 = plot(xSREmpirical,ySREmpirical)
    ax.set_title('Success Rate (100 iterations)', fontsize=16)
    ylabel('Time Taken (s)', fontsize=18)
    xlabel('Number of Measurements', fontsize=15)
    axis([5, 5000, 0, 3330])

    ax = subplot(132)
    p1 = plot(xPI, yPI)
    ax.set_title('Perceived Information', fontsize=15)
    xlabel('Number of Measurements', fontsize=15)

    ax = subplot(133)
    p1 = plot(xSNR, ySNR)
    ax.set_title('Signal to Noise Ratio', fontsize=15)
    legend(loc='lower right')
    xlabel('Number of Measurements', fontsize=15)

    show()

