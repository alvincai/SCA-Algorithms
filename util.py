#! /usr/bin/python

def hex2str(a):
    """Return a string of an 8-bit hex value. For example:
        0xf1  -> 'f1'
        0x2   -> '02'
        0x234 -> '34'
    """
    a &= 0xff
    a = int(a)
    s = ""
    if (a < 0x10):
        s = "0"
    s += hex(a)
    s = s.replace("0x", "")
    s = s.replace("L", "")
    return s

def HW(num):
    """Returns hamming weight of num"""
    return bin(num).count("1")

def HD(num1, num2):
    """Returns hamming distance of num1 and num2"""
    return HW(num1^num2)


