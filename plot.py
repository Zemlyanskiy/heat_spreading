# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import re
import sys

PATH=sys.argv[1]

#--------------------------------------READING FILE----------------------------------------
file=open(PATH).readlines()
matrix=[]
for line in file:
    if re.match(r"t=",line):
        t=float(re.search(r"[\d.]+",line).group())
    elif re.match(r"T=",line):
        T=float(re.search(r"[\d.]+",line).group())
    elif re.match(r"deltaT=",line):
        deltaT=float(re.search(r"[\d.]+",line).group())
    elif re.match(r"XMin=",line):
        XMin=float(re.search(r"[-\d.]+",line).group())
    elif re.match(r"XMax=",line):
        XMax=float(re.search(r"[-\d.]+",line).group())
    elif re.match(r"Lx=",line):
        Lx=int(re.search(r"[\d.]+",line).group())
    elif re.match(r"Sigma=",line):
        Sigma=float(re.search(r"[\d.]+",line).group())
    elif re.match(r"deltaOut=",line):
        deltaOut=float(re.search(r"[\d.]+",line).group())
    else:
        matrix.append(re.findall(r"[\d.]+",line))

step = (XMax-XMin)/Lx
i=XMin
xArr=[]
globals()['iterator']=1
#---------------------------------START DRAWING PLOT------------------------------------------
freqs = np.arange(2, 20, 3)

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)
t = np.arange(0.0, 1.0, 0.001)
s = np.sin(2*np.pi*freqs[0]*t)

xcoords=np.arange(XMin,XMax,(XMax-XMin)/Lx)
l, = plt.plot(xcoords,matrix[0], lw=2)
plt.xlabel('OX')
plt.ylabel('OY')
ax.set_title('0')
class Index(object):
    ind = 0

    def next(self, event):
        if self.ind < len(matrix)-1:
            self.ind += 1
            i = self.ind % len(matrix)
            ydata = matrix[i]
            l.set_ydata(ydata)
            ax.set_title(self.ind)
            plt.draw()

    def prev(self, event):
        if self.ind>0:
            self.ind -= 1
            i = self.ind % len(matrix)
            ydata = matrix[i]
            l.set_ydata(ydata)
            ax.set_title(self.ind)
            plt.draw()

callback = Index()
axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)

plt.show()
exit(0)