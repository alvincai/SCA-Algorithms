#! /usr/bin/python

from numpy import sin, linspace, pi
from pylab import plot, show, title, xlabel, ylabel, subplot
from scipy import fft, arange
import LoadTraces

def plotSpectrum(y,Fs):
    """
    Plots a Single-Sided Amplitude Spectrum of y(t)
    """
    n = len(y) # length of the signal
    k = arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range

    Y = fft(y)/n # fft computing and normalization
    Y = Y[range(n/2)]

    plot(frq,abs(Y),'r') # plotting the spectrum
    xlabel('Freq (MHz)')
    ylabel('|Y(freq)|')

if __name__ == "__main__":

    l = LoadTraces.LoadTraces('traces/0004.trs', numTraces = 2)
    plotSpectrum(l.traces[1], 100.0)
    show()
