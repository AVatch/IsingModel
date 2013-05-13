import numpy as np
import matplotlib.pyplot as plt

fileName = 'Ising2D_Data.dat'
#set up fig
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2)

#import data
T, M_avg, Mabs_avg, Msq_avg, x, xabs, E_avg, Esq_avg, C, Ul = np.loadtxt(fileName, unpack=True)

#plot M_avg against T
ax1.set_title('M_avg')
ax1.plot(T,M_avg,'k')
ax1.plot(T,M_avg,'k-')

#Plot E_avg
ax2.set_title('E_avg')
ax2.plot(T,E_avg,'b')
ax2.plot(T,E_avg,'b-')

#Plot Heat Capacity
ax3.set_title('C')
ax3.plot(T,C,'g')
ax3.plot(T,C,'g-')

#Plot susceptibility
ax4.set_title('x')
ax4.plot(T,x,'r')
ax4.plot(T,x,'r-')

plt.tight_layout()
plt.show()

 
