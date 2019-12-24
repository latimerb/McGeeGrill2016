import numpy as np
import h5py
import pdb

vars_file = './output/cell_vars.h5'

f = h5py.File(vars_file,'r')
plt.plot(f['report']['LUT']['data'][:,1])
plt.show()

