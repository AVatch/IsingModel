import numpy as np
import matplotlib.pyplot as plt

#import data
T, M_avg, Mabs_avg, Msq_avg, x, xabs, E_avg, Esq_avg, C, Ul = np.loadtxt('Ising2D_Data_T4_1.dat', unpack=True)

plt.plot(T,E_avg,'r-',label='E_avg')
plt.plot(T,M_avg,'b--',label='M_avg')
plt.grid()
plt.legend()
plt.xlabel('Temp (K)')
plt.title('Phase transition around 1.8K')

plt.show()
