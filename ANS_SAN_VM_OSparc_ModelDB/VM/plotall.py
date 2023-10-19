import sys, pprint
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from io import StringIO   # StringIO behaves like a file object



ap1 = np.loadtxt('ap_icnsach.txt');
#column 1: time
#column 2: voltage
                
cai1 = np.loadtxt('ca_icnsach.txt');
#column 1: time
#column 2: ca concentration

#figure;
fig3, axs3 = plt.subplots(2, 1, figsize=(12.8, 9.6) )
fig3.subplots_adjust(wspace=0.5)

axs3[0].plot( ap1[:,0]*0.001, ap1[:,1], 'm' )
axs3[0].set_ylim( [-100, 50] )

axs3[0].set_xlabel('Time (sec)', fontsize=20)
axs3[0].set_ylabel('Vm (mV)', fontsize=20)


axs3[1].plot( cai1[:,0]*0.001, cai1[:,1]*1000, 'm' )
axs3[1].set_ylim( [0, 1] )

axs3[1].set_xlabel('Time (sec)', fontsize=20)
axs3[1].set_ylabel(r'Cai [$\mu$M]', fontsize=20)

fig3.savefig( 'figpython.png' )

#
#plt.show()
