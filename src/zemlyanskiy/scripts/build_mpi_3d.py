import sys
import os
import shutil

arg_1 = sys.argv[1]#code filename
arg_2 = sys.argv[2]#txt filename
arg_3 = sys.argv[3]#kind of method
arg_4 = sys.argv[4]#choice of method
arg_5 = sys.argv[5]#number of threads
arg_6 = sys.argv[6]#number of processes

os.system(r'mpicc' +' '+ arg_1 +' '+ '-o build -m64 -O2 -ftree-vectorize -flto -march=native -funroll-loops -lm -fopenmp')
os.system(r'mpirun -np' + ' ' + arg_6 +' ' + './build'+' '+
            arg_2 +' ' + arg_3 +' '+ arg_4 + ' ' + arg_5)
os.system(r'python3 ..//..//heat_map.py'+' '+arg_2)
    