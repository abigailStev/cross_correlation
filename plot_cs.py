import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

font_prop = font_manager.FontProperties(size=18)
table = np.loadtxt('cs_avg_adjusted.dat')
table = np.loadtxt('cs_avg.dat')

freq = np.asarray(table[:,0], dtype=np.float)
cs_real = np.asarray(table[:,1], dtype=np.float)

pos_freq = freq[np.where(freq >= 0)]
pos_cs_real = cs_real[np.where(freq >= 0)]

fig,ax = plt.subplots(1,1,figsize=(10,8))
ax.plot(pos_freq, pos_cs_real, lw=2)
ax.hlines(0.0, 4, 7, linestyle='dashed', lw=1.5, color='gray')
ax.vlines(5.42176, -5e7, 5e7, linestyle='dotted', lw=1.5, color='gray')
ax.set_xlim(4,7)
ax.set_ylim(-5e7, 5e7)
ax.set_xlabel("Frequency (Hz)", fontproperties=font_prop)
ax.set_ylabel("Amplitude of real component of cross spectrum", fontproperties=font_prop)
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
plt.title("Cross Spectrum", fontproperties=font_prop)

plt.show()
plt.savefig('cs_avg.png', dpi=200)
plt.close()
