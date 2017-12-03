# -*- coding: utf8 -*-

import re
import matplotlib.pyplot as plt
file=open("input1.txt").readlines()
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
while i<XMax:
    xArr.append(i)
    i+=step
result=matrix[0]
del matrix[0]
line=plt.plot(xArr,result)
print("end step\n")
for result in matrix:
    line[0].set_ydata(result)
    print("end step\n")
exit(0)