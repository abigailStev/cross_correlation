import numpy as np
import matplotlib.pyplot as plt

table = np.loadtxt('cs_avg.dat')

cs = np.asarray(table[:,0])
fig = plt.subplots()
bins = np.arange(len(cs))
plt.plot(bins, cs)
plt.show()

