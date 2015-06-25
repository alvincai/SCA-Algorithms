#! /usr/bin/python
# Convert Inspector traceset to numpy array.
# Author: Alvin

import numpy as np
import matplotlib.pyplot as plt
import random
import gc
import random
import struct
import binascii

class LoadTraces():
    """Converts Inspector traceset to numpy array
    This is required for other classes (SNR, DES_CPA, Success rate etc) to function
    """

    def __init__(self, fname, numTraces=0, samplePoint =None, offset = 0):
        """
        Loads the inspector traceset

        #Inputs:
        fname - file name of the inspector traceset
        numTraces - Optional, will load all traces if unspecified, otherwise just the specified numTraces in the file after the offset
        samplePoint - Optional, will load all samples in the traceset if unspecified. Otherwise just the specified sample point
        Offset - See description for numTraces


        """
        self.offset = 0
        self.directory = ""
        self.numTraces = 0
        self.numTracesUser = numTraces
        self.numSamples = 0
        self.floatingPoint = True    # if false, data encoded as Integer
        self.bytesPerSample = 1
        self.sampleType = None
        self.cryptoDataLength = 0
        self.blockSize = 0          # Block size of the plaintext/ ciphertext. Should be 1/2 of cryptoDataLength
        self.titleSpace = 0
        self.trsData = bytearray() # Data of the trs file
        self.log = []
        self.f = None           # file handler
        self.endOfHeader = False
        self.inputtext = None
        self.outputtext = None
        self.traces = None
        self.samplePoint = samplePoint

        self.headerCoding = { '\x41' : self.setNumTraces,    # Number of Traces
                     '\x42' : self.setNumSamples,   # Number of Samples per trace
                     '\x43' : self.setSampleCoding, # Sample Coding e.g. integer, floating point, length in bytes
                     '\x44' : self.setCryptoDataLength,  # length of cryptographic data included in trace
                     '\x45' : self.setTitleSpace,        # Title space reserved per trace
                     '\x46' : self.ignore,          # Global trace title
                     '\x47' : self.ignore,          # Description
                     '\x48' : self.ignore,          # Offset in X-axis for trace representation
                     '\x49' : self.ignore,          # Label of X-axis
                     '\x4A' : self.ignore,          # Label of Y-axis
                     '\x4B' : self.ignore,          # Scale value for X-axis
                     '\x4C' : self.yscale,          # Scale value for Y-axis
                     '\x4D' : self.ignore,          # Trace offset for displaying trace number
                     '\x4E' : self.ignore,          # Log Scale
                     '\x4F' : self.ignore,          # Reserved
                     '\x50' : self.ignore,          # Reserved
                     '\x51' : self.ignore,          # Reserved
                     '\x52' : self.ignore,          # Reserved
                     '\x53' : self.ignore,          # Reserved
                     '\x54' : self.ignore,          # Reserved
                     '\x55' : self.ignore,          # Range of the scope
                     '\x56' : self.ignore,          # Coupling of the scope
                     '\x57' : self.ignore,          # Offset of the scope
                     '\x58' : self.ignore,          # Input impedance of the scope
                     '\x59' : self.ignore,          # Device ID of the scope
                     '\x5A' : self.ignore,          # Type of filter used
                     '\x5B' : self.ignore,          # Frequency of the filter used
                     '\x5C' : self.ignore,          # Range of the filter used
                     '\x5F' : self.endHeader,       # Trace block marker that indicates end of the header
                     '\x60' : self.ignore,          # Ext clock used
                     '\x61' : self.ignore,          # Ext clock threshold
                     '\x62' : self.ignore,          # Ext clock multiplier
                     '\x63' : self.ignore,          # Ext clock phase shift
                     '\x64' : self.ignore,          # Ext clock resampler mask
                     '\x65' : self.ignore,          # Ext clock resampler enabled
                     '\x66' : self.ignore,          # Ext clock freq
                     '\x67' : self.ignore,          # Ext clock time base
        }

        self.load(fname)



    def load(self, filename):
        self.f = open(filename, "rb")
        try:
            byte = self.f.read(1)

            # Read Inspector Header of TLV format
            while (self.endOfHeader == False):
                self.headerCoding[byte]()
                byte = self.f.read(1)
            if (self.numTracesUser != 0) and (self.numTracesUser <self.numTraces):
                self.numTraces = self.numTracesUser

            # Read trace data
            self.readINSdata()
        finally:
            self.f.close()


    def readINSdata(self):
        """
        Read binary data that consists of -
            numTraces times of:
                $titleSpace bytes of title space (ignored)
                $cryptoDateLength bytes of plaintext and ciphertext data (saved to inputdata.npy and outputdata.npy)
                $numSamples of int/float samples of bytesPerSample bytes (saved to trace.npy)
        """
        if (self.cryptoDataLength % 2 == 0):
            self.blockSize = self.cryptoDataLength / 2
        else:
            print "Error: length of crypto data is odd"
            return

        self.inputtext = np.zeros((self.numTraces, self.blockSize), dtype = np.uint8)
        self.outputtext = np.zeros((self.numTraces, self.blockSize), dtype = np.uint8)

        if (self.floatingPoint and  (self.bytesPerSample == 1)):
            print "Error: unlikely that inspector float is 8 bit..."
        elif (self.floatingPoint and  (self.bytesPerSample == 2)):
            self.sampleType = np.float16
        elif (self.floatingPoint and  (self.bytesPerSample == 4)):
            self.sampleType = np.float32
        elif ((self.floatingPoint == False) and  (self.bytesPerSample == 1)):
            self.sampleType = np.byte
        elif ((self.floatingPoint == False) and  (self.bytesPerSample == 2)):
            self.sampleType = np.short
        elif ((self.floatingPoint == False) and  (self.bytesPerSample == 4)):
            self.sampleType = np.int32

        if (self.samplePoint == None):
            self.traces = np.zeros((self.numTraces, self.numSamples), dtype = self.sampleType)
        else:
            self.traces = np.zeros((self.numTraces, 1), dtype = self.sampleType)
            tmpTrace = np.zeros(self.numSamples, dtype = self.sampleType)

        # Discard the first few traces until the specified offset
        for i in range(self.offset):
            self.f.read(self.titleSpace)
            self.f.read(self.blockSize)
            self.f.read(self.blockSize)
            self.f.read(self.numSamples * self.bytesPerSample)


        for i in range(self.numTraces):
            self.f.read(self.titleSpace)    #discard
            ba = self.f.read(self.blockSize)
            self.inputtext[i] = np.fromstring(ba, dtype = np.uint8)

            ba = self.f.read(self.blockSize)
            self.outputtext[i] = np.fromstring(ba, dtype = np.uint8)

            ba = self.f.read(self.numSamples * self.bytesPerSample)
            if (self.samplePoint == None):
                self.traces[i] = np.fromstring(ba, dtype = self.sampleType)
            else:
                tmpTrace = np.fromstring(ba, dtype = self.sampleType)
                self.traces[i] = tmpTrace[self.samplePoint]

        if (self.samplePoint != None):
            self.numSamples = 1



    def readTLV(self):
        """ Read Length and Value of header object, then returns value as bytearray"""
        byte = self.f.read(1)
        length = int(byte.encode('hex'), 16)
        ba = self.f.read(length)
        return ba


    def setNumTraces(self):
        ba = self.readTLV()
        for i in range(4):
            value = int(ba[i].encode('hex'), 16)
            self.numTraces |= ((value & 0xff) << (i * 8))
        print ("Number of Traces is %d" % self.numTraces)


    def setNumSamples(self):
        ba = self.readTLV()
        for i in range(4):
            value = int(ba[i].encode('hex'), 16)
            self.numSamples |= ((value & 0xff) << (i * 8))
        print ("Number of Samples per Trace is %d" % self.numSamples)


    def yscale(self):
        ba = self.readTLV()
        #print binascii.hexlify(ba)
        self.yscale = struct.unpack('f', ba)
        self.yscale = self.yscale[0]
        print ("yscale is %f" % self.yscale)


    def setSampleCoding(self):
        ba = self.readTLV()
        value = int(ba[0].encode('hex'), 16)
        if (value & 0x10):
            self.floatingPoint = True
            print "Each sample is stored as floating point, of length ",
        else:
            self.floatingPoint = False
            print "Each sample is stored as Integer, of length ",

        self.bytesPerSample = value & 0x0f
        print ("%d bytes" % self.bytesPerSample)


    def setCryptoDataLength(self):
        ba = self.readTLV()
        for i in range(2):
            value = int(ba[i].encode('hex'), 16)
            self.cryptoDataLength |= ((value & 0xff) << (i * 8))
        print "Crypto Data is of length %d bytes" % self.cryptoDataLength


    def setTitleSpace(self):
        ba = self.readTLV()
        self.titleSpace = int(ba[0].encode('hex'), 16)


    def ignore(self):
        self.readTLV()

    def endHeader(self):
        self.endOfHeader = True


    def addNoise(self, noiseRange, samplePoint):
        for i in range(self.numTraces):
            noise = random.random() * noiseRange
            self.traces[i][samplePoint] += noise


    def reduceTrace(self, samplePoint):
        """Delete all but the selected samplepoint in traceset file"""

        self.numSamples = 1
        tracesNew = np.zeros((self.numTraces, self.numSamples), dtype = self.sampleType)
        for i in range(self.numTraces):
            tracesNew[i][0] = self.traces[i][samplePoint]
        self.traces = tracesNew
        gc.collect()

    def shuffle(self):
        """Randomly shuffle the traces"""
        tracesNew = np.zeros((self.numTraces, self.numSamples), dtype = self.sampleType)
        inputtextNew = np.zeros((self.numTraces, self.blockSize), dtype = np.uint8)
        outputtextNew = np.zeros((self.numTraces, self.blockSize), dtype = np.uint8)

        for i in range (self.numTraces):
            rand = random.randint(0, self.numTraces - i -1)
            tracesNew[i] = self.traces[rand]
            inputtextNew[i] = self.inputtext[rand]
            outputtextNew[i] = self.outputtext[rand]
            self.traces = np.delete(self.traces, rand, axis=0)
            self.inputtext = np.delete(self.inputtext, rand, axis=0)
            self.outputtext = np.delete(self.outputtext, rand, axis=0)

        self.traces = tracesNew
        self.inputtext = inputtextNew
        self.outputtext = outputtextNew
        gc.collect()

if __name__ == "__main__":
    l = LoadTraces('traces/2001.trs', numTraces = 2)

    plt.plot(l.traces[0])
    plt.show()

