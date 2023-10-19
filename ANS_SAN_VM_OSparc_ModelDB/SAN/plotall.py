import sys, pprint
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from io import StringIO   # StringIO behaves like a file object



readallach = np.loadtxt( 'alloutputs_icnsach.txt' )
tach = readallach[:,0]/1e3 # time
vmach = readallach[:,15] # voltage


readICach = np.loadtxt('outputs_icnsach.txt')
tcach = readICach[:,0]/1e3  # time

cch = readICach[:,5] # Ach
icns = readICach[:,6] # icns



mycolor = 'b';


#figure;
fig3, axs3 = plt.subplots(2, 1, figsize=(12.8, 9.6) )
fig3.subplots_adjust(wspace=0.5)
#subplot(1,3,1)
axs3[0].plot( tcach, cch, 'g' )
axs3[0].plot( tcach, icns, 'b' )
axs3[0].set_ylim( [0, 150] )

# calculate heart rates
loc, properties = find_peaks( vmach , prominence=0.01, width=20 )
HR = 60./np.diff(loc) * 1000

# plot heart rate
axs3[1].plot( loc[1:]*0.001, HR, 'm')
axs3[1].set_xlabel('Time (sec)', fontsize=20)
axs3[1].set_ylabel('AP firing rate [bpm]', fontsize=20)

fig3.savefig( 'figpython.png' )

#
#plt.show()
