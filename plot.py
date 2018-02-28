# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
import re
import sys

def check_abs_error(diff, abs_err):
    flag = 1
    for element in diff:
        if element > abs_err:
            flag=0
    if flag==1 :
        fig.suptitle('OK', fontsize=15)
    if flag==0 :
        fig.suptitle('ERROR', fontsize=15)
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

step = (file1.XMax-file1.XMin)/file1.Lx
i=file1.XMin
xArr=[]
globals()['iterator']=1
abs_err=0.1
#---------------------------------START DRAWING PLOT------------------------------------------
freqs = np.arange(2, 20, 3)
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)
t = np.arange(0.0, 1.0, 0.001)
s = np.sin(2*np.pi*freqs[0]*t)
xcoords=np.arange(file1.XMin,file1.XMax,(file1.XMax-file1.XMin)/file1.Lx)
l,=plt.plot(xcoords,file1.matrix[0])
plt.xlabel('OX')
plt.ylabel('OY')
ax.set_title('0')

#FUNCTIONALUTY FOR 2 FILE
if len(sys.argv)==3:
    PATH2 = sys.argv[2]
    file2 = FILE_FORMAT(PATH2)
    l2, = plt.plot(xcoords, file2.matrix[0])
    difference=[];
    buf_list=[];
    for el_file1,el_file2 in zip(file1.matrix,file2.matrix):
        for num_el1,num_el2 in zip(el_file1,el_file2):
            buf_list.append(abs(num_el1-num_el2))
        difference.append(buf_list)
        buf_list=[]
    l3, = plt.plot(xcoords, difference[0])
    plt.xlabel('first is green, second is orange, difference is green')
    check_abs_error(difference[0],abs_err)





#BUTTONS FUNCTIONALITY
class Index(object):
    ind = 0
    def next(self, event):
        if self.ind < len(file1.matrix)-1:
            self.ind += 1
            i = self.ind % len(file1.matrix)
            data=file1.matrix[i]
            ydata = data
            if len(sys.argv) == 3:
                data2 = file2.matrix[i]
                ydata2=data2
                l2.set_ydata(ydata2)
                diff_data2=difference[i]
                ydata_diff=diff_data2
                l3.set_ydata(ydata_diff)
                check_abs_error(difference[i], abs_err)
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
            if len(sys.argv) == 3:
                data2 = file2.matrix[i]
                ydata2=data2
                l2.set_ydata(ydata2)
                diff_data2 = difference[i]
                ydata_diff = diff_data2
                l3.set_ydata(ydata_diff)
                check_abs_error(difference[i], abs_err)
            ax.set_title(self.ind)
            plt.draw()
    def submit(self, text):
        abs_err=float(text)
        check_abs_error(difference[self.ind], abs_err)
callback = Index()
axprev = plt.axes([0.7, 0.01, 0.06, 0.075])
axnext = plt.axes([0.81, 0.01, 0.06, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next,)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)

#---------------------------------DRAWING TEXT FIELDS------------------------------------------
initial_text = "0.1"
axbox = plt.axes([0.40, 0.01, 0.06, 0.055])
text_box = TextBox(axbox, 'Abs', initial=initial_text)
text_box.on_submit(callback.submit)

plt.show()
exit(0)