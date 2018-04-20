# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
import re
import sys

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
# --------------------------------------READING FILE----------------------------------------

PATH = sys.argv[1]
file1=FILE_FORMAT(PATH)
#---------------------------------START DRAWING PLOT------------------------------------------

freqs = np.arange(2, 20, 3)
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)
xcoords=np.arange(file1.XMin,file1.XMax,(file1.XMax-file1.XMin)/file1.Lx)
l,=plt.plot(xcoords,file1.matrix[0])
plt.xlabel('OX')
plt.ylabel('OY')
ax.set_title('0')


#BUTTONS FUNCTIONALITY
class Index(object):
    ind = 0
    def next(self, event):
        if self.ind < len(file1.matrix)-1:
            self.ind += 1
            i = self.ind % len(file1.matrix)
            data=file1.matrix[i]
            ydata = data
            l.set_ydata(ydata)
            ax.set_title(self.ind)
            plt.draw()
    def prev(self, event):
        if self.ind>0:
            self.ind -= 1
            i = self.ind % len(file1.matrix)
            data=file1.matrix[i]
            ydata = data
            l.set_ydata(ydata)
            ax.set_title(self.ind)
            plt.draw()
callback = Index()
axprev = plt.axes([0.5, 0.01, 0.2, 0.075])
axnext = plt.axes([0.71, 0.01, 0.2, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next,)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)

plt.show()
exit(0)