import numpy as np
import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Usage: python plotter.py `file.out`")
    exit(0)

x = np.fromfile(sys.argv[1], np.float64)
n = int(len(x) **.5)
x = x.reshape((n,n))
plt.imshow(x)
plt.show()
