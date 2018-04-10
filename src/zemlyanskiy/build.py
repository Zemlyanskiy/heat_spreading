import sys
import os
import shutil

arg_1 = sys.argv[1]
arg_2 = sys.argv[2]
arg_3 = sys.argv[3]
arg_4 = sys.argv[4]

os.system(r'gcc' +' '+ arg_1 +' '+ '-o build -lm -fopenmp')
os.system(r'./build'+' '+arg_2 +' ' + arg_3 +' ' + arg_4)
if len(sys.argv)==6:
    arg_5 = sys.argv[5]
    os.system(r'python3 ..//..//plot.py'+' '+arg_2 + ' '+ '..//..//output/zemlyansky/'+arg_5)
else:
    os.system(r'python3 ..//..//plot.py'+' '+arg_2)
if arg_3 == 'e':
    shutil.copyfile(arg_2, r'..//..//output/zemlyansky/eiler_output.txt')
elif arg_3 == 'r':
    shutil.copyfile(arg_2, r'..//..//output/zemlyansky/runge_output.txt')
