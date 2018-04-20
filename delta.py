# -*- coding: utf8 -*-
import re
import sys
import numpy as np


class FILE_FORMAT:
    t=None
    T=None
    deltaT=None
    XMin=None
    XMax=None
    Lx=None
    Sigma=None
    deltaOut=None
    matrix = None
    def __init__(self, PATH):
        file = open(PATH).readlines()
        self.matrix = []
        for line in file:
            if re.match(r"t=", line):
                self.t = float(re.search(r"[\d.]+", line).group())
            elif re.match(r"T=", line):
                self.T = float(re.search(r"[\d.]+", line).group())
            elif re.match(r"deltaT=", line):
                self.deltaT = float(re.search(r"[\d.]+", line).group())
            elif re.match(r"XMin=", line):
                self.XMin = float(re.search(r"[-\d.]+", line).group())
            elif re.match(r"XMax=", line):
                self.XMax = float(re.search(r"[-\d.]+", line).group())
            elif re.match(r"Lx=", line):
                self.Lx = int(re.search(r"[\d.]+", line).group())
            elif re.match(r"Sigma=", line):
                self.Sigma = float(re.search(r"[\d.]+", line).group())
            elif re.match(r"deltaOut=", line):
                self.deltaOut = float(re.search(r"[\d.]+", line).group())
            else:
                tmp = re.findall(r"[\d.]+", line)
                float_list = []
                for el in tmp:
                    float_list.append(float(el))
                self.matrix.append(float_list)
                float_list = []

ABS=0.001
PATH = sys.argv[1]
file1=FILE_FORMAT(PATH)

PATH2 = sys.argv[2]
file2=FILE_FORMAT(PATH2)

ABS= float(sys.argv[3])

FLAG=0;

if ((file1.t == file2.t)and(file1.T == file2.T)and(file1.deltaT == file2.deltaT)and(file1.XMax == file2.XMax)and(file1.XMin == file2.XMin)and(file1.XMax == file2.XMax)and(file1.Lx == file2.Lx)and(file1.deltaOut == file2.deltaOut)and(file1.Sigma == file2.Sigma)):
    for line_file1,line_file2 in zip(file1.matrix, file2.matrix):
        for el_file1, el_file2 in zip (line_file1,line_file2):
            if (abs(el_file1-el_file2)>ABS):
                FLAG+=1
else:
    print("\nNON COMPATIBLE FILES\n")
    exit(2)

print("\nNUMBER OF ERROR " + str(FLAG)+"\n")
exit(0)

