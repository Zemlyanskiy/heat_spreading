# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
PATH2 = sys.argv[2]
file1=FILE_FORMAT(PATH)
file2 = FILE_FORMAT(PATH2)
difference=[]
relative_difference=[]
buf_list=[]
for el_file1,el_file2 in zip(file1.matrix,file2.matrix):
    for num_el1,num_el2 in zip(el_file1,el_file2):
        buf_list.append(abs(num_el1-num_el2))
    difference.append(buf_list)
    buf_list=[]

#for line_file1,line_file2,line_diff in zip(file1.matrix,file2.matrix,difference):
#   for el_file1,el_file2,el_dif in zip(line_file1,line_file2,line_diff):

#---------------------------------START DRAWING PLOT------------------------------------------
xcoords=np.arange(file1.XMin,file1.XMax,(file1.XMax-file1.XMin)/file1.Lx)

# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 2)

fig = plt.figure()

ax1 = fig.add_subplot(gs[0, 0]) # row 0, col 0
ax1.plot(xcoords, file1.matrix[0], color="blue")

ax2 = fig.add_subplot(gs[0, 1]) # row 0, col 1
ax2.plot(xcoords, file2.matrix[0], color="red")

axabs = fig.add_subplot(gs[1,:]) # row 1, span all columns
axabs.plot(xcoords, difference[0], color="brown")

ax1.set_ylim(0, 1)
ax2.set_ylim(0, 1)
axabs.set_ylim(0, 1)
ax1.set_title("First File", fontsize=7)
ax2.set_title("Second File", fontsize=7)
axabs.set_title("Absolete", fontsize=7)
#BUTTONS FUNCTIONALITY
class Index(object):
    abs_err = 0.1
    ind = 0
    def check_abs_error(self,diff):
        flag = 1
        for element in diff:
            if element > self.abs_err:
                flag = 0
        if flag == 1:
            fig.suptitle(str(self.ind) + ' OK', fontsize=15)
        if flag == 0:
            fig.suptitle(str(self.ind) + ' ERROR', fontsize=15)
    def PreSet(self):
        ax1.clear()
        ax2.clear()
        axabs.clear()
        ax1.set_ylim(0, 1)
        ax2.set_ylim(0, 1)
        axabs.set_ylim(0, 1)
        ax1.set_title("First File", fontsize=7)
        ax2.set_title("Second File", fontsize=7)
        axabs.set_title("Absolete", fontsize=7)
    def next(self, event):
        if self.ind < len(file1.matrix)-1:
            self.PreSet()
            self.ind += 1
            i = self.ind % len(file1.matrix)

            ax1.plot(xcoords, file1.matrix[i], color="blue")
            ax2.plot(xcoords, file2.matrix[i], color="red")
            axabs.plot(xcoords, difference[i], color="orange")


            self.check_abs_error(difference[i])

            plt.draw()
    def prev(self, event):
        if self.ind > 0:
            self.PreSet()
            self.ind -= 1
            i = self.ind % len(file1.matrix)

            ax1.plot(xcoords, file1.matrix[i], color="blue")
            ax2.plot(xcoords, file2.matrix[i], color="red")
            axabs.plot(xcoords, difference[i], color="orange")
            self.check_abs_error(difference[i])

            plt.draw()

    def submit(self, text):
        self.abs_err=float(text)
        self.check_abs_error(difference[self.ind])


callback = Index()
axprev = plt.axes([0.5, 0.01, 0.2, 0.045])
axnext = plt.axes([0.71, 0.01, 0.2, 0.045])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next,)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)


initial_text = "0.1"
axbox = plt.axes([0.20, 0.01, 0.20, 0.045])
text_box = TextBox(axbox, 'Abs', initial=initial_text)
text_box.on_submit(callback.submit)

plt.show()
exit(0)