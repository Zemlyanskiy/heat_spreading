# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button, TextBox

import re
import sys

class FileFormat:
    t=None
    T=None
    deltaT=None
    XMin=None
    XMax=None
    YMin=None
    YMax=None
    ZMin=None
    ZMax=None
    Lx=None
    Ly=None
    Lz=None
    Sigma=None
    deltaOut=None
    data = []
    def __init__(self, PATH):
        file = open(PATH).readlines()
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
            elif re.match(r"YMin=", line):
                self.YMin = float(re.search(r"[-\d.]+", line).group())
            elif re.match(r"YMax=", line):
                self.YMax = float(re.search(r"[-\d.]+", line).group())
            elif re.match(r"ZMin=", line):
                self.ZMin = float(re.search(r"[-\d.]+", line).group())
            elif re.match(r"ZMax=", line):
                self.ZMax = float(re.search(r"[-\d.]+", line).group())
            elif re.match(r"Lx=", line):
                self.Lx = int(re.search(r"[\d.]+", line).group())
            elif re.match(r"Ly=", line):
                self.Ly = int(re.search(r"[\d.]+", line).group())
            elif re.match(r"Lz=", line):
                self.Lz = int(re.search(r"[\d.]+", line).group())
            elif re.match(r"Sigma=", line):
                self.Sigma = float(re.search(r"[\d.]+", line).group())
            elif re.match(r"deltaOut=", line):
                self.deltaOut = float(re.search(r"[\d.]+", line).group())
            else:
                data_snapshot =  np.zeros((self.Lx, self.Ly, self.Lz))
                tmp = re.findall(r"[\d.]+", line)

                for i in range(0,self.Lx):
                    for j in range(0,self.Ly):
                        for k in range(0,self.Lz):
                            #if float(tmp[i*self.Lx*self.Ly+j*self.Ly+k]) != 0:
                            #    print(i,j,k)
                            data_snapshot[i][j][k] = float(tmp[i*self.Lx*self.Ly+j*self.Ly+k])
                self.data.append(data_snapshot)
# --------------------------------------READING FILE----------------------------------------

file= FileFormat(sys.argv[1])

section_number = int(file.Lz/2)
# generate 2 2d grids for the x & y bounds
y, x = np.meshgrid(np.linspace(file.XMin, file.XMax, file.Lx), np.linspace(file.YMin, file.YMax, file.Ly))
z = np.zeros((file.Lx,file.Ly))
for i in range(0, file.Lx):
    for j in range(0, file.Ly):
            z[i][j] = file.data[0][i][j][section_number]

# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
z = z[:-1, :-1]
z_min, z_max = -np.abs(z).max(), np.abs(z).max()
fig, ax = plt.subplots()

c = ax.pcolormesh(x, y, z, cmap='coolwarm', vmin=z_min, vmax=z_max)
ax.set_title("Section: %d, Time: %f" % (section_number, 0))
# set the limits of the plot to the limits of the data
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax=ax)

#BUTTONS FUNCTIONALITY
class Index(object):
    index = int(file.Lz/2)
    current_time = 0
    def PreSet(self):
        ax.clear()
        ax.set_title("Section: %d, Time: %f" % (self.index, self.current_time))
        ax.axis([x.min(), x.max(), y.min(), y.max()])

    def ReDraw(self):
        for i in range(0, file.Lx - 1):
            for j in range(0, file.Ly - 1):
                z[i][j] = file.data[self.current_time][i][j][self.index]

        c = ax.pcolormesh(x, y, z, cmap='coolwarm', vmin=z_min, vmax=z_max)
        plt.draw()

    def nextSection(self, event):
        if self.index < len(file.data[0])-1:
            self.PreSet()
            self.index += 1
            self.ReDraw()

    def prevSection(self, event):
        if self.index > 0:
            self.PreSet()
            self.index -= 1
            self.ReDraw()

    def nextTime(self, event):
        if self.current_time < len(file.data)-1:
            self.PreSet()
            self.current_time += 1
            self.ReDraw()

    def prevTime(self, event):
        if self.current_time > 0:
            self.PreSet()
            self.current_time -= 1
            self.ReDraw()

callback = Index()

axprevtime = plt.axes([0.1, 0.01, 0.2, 0.045])
bprevtime = Button(axprevtime, 'Previous Time')
bprevtime.on_clicked(callback.prevTime)

axnexttime = plt.axes([0.3, 0.01, 0.2, 0.045])
bnexttime = Button(axnexttime, 'Next Time')
bnexttime.on_clicked(callback.nextTime)

axprevsec = plt.axes([0.5, 0.01, 0.2, 0.045])
bprevsec = Button(axprevsec, 'Previous Section')
bprevsec.on_clicked(callback.prevSection)

axnextsec = plt.axes([0.71, 0.01, 0.2, 0.045])
bnextsec = Button(axnextsec, 'Next Section')
bnextsec.on_clicked(callback.nextSection)

plt.show()