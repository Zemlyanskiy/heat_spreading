import sys
import os
import shutil

arg_1 = sys.argv[1]#code filename
os.system(r'gcc' +' '+ arg_1 +' '+ '-g -lm -fopenmp')