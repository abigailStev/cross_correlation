import numpy as np
import matplotlib.pyplot as plt

table = np.loadtxt('cs_avg.dat')

freq = np.asarray(table[:,0], type=np.complex128)
cs = np.asarray(table[:,1], type=np.complex128)
fig,ax = plt.subplots(1,1,figsize=(6,6))
ax.plot(freq, cs)
ax.set_xlim(0,20)
plt.show()

