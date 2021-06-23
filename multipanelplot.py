import matplotlib.pyplot as plt
import numpy as np

x,y = np.loadtxt('ruth_1.75_data.data', delimiter=',', unpack=True)
x2,y2 = np.loadtxt('ruth_3.7_data.data', delimiter=',', unpack=True)



fig, axs = plt.subplots(2,2)

axs[0, 0].plot(x, y, 'tab:blue')
#axs[0, 0].set_title('Axis [0, 0]')
axs[0,0].set_ylabel(r'log $\rho_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
#axs[0,0].legend('z = 2', handlelength=2)
axs[0,0].text(-1.6,8.0,  'z = 1.75', fontsize=12)
m,b = np.polyfit(x,y,1)
axs[0,0].plot(x,m*x+b,color='black')


axs[0, 1].plot(x2, y2, 'tab:orange')
#axs[0, 1].set_title('Axis [0, 1]')
axs[0,1].text(-1.6,7.7, 'z = 3.79', fontsize=12)
m,b = np.polyfit(x2,y2,1)
axs[0,1].plot(x2,m*x2+b,color='black')

"""
axs[1, 0].plot(x3, y3, 'tab:green')
#axs[1, 0].set_title('Axis [1, 0]')
axs[1,0].set_ylabel('y-label')
axs[1,0].set_xlabel('x-label')
axs[1,0].text(225,600, 'z = 4', fontsize=12)

axs[1, 1].plot(x4, y4, 'tab:red')
#axs[1, 1].set_title('Axis [1, 1]')
axs[1,1].set_xlabel('x-label')
axs[1,1].text(150,450, 'z = 5', fontsize=12)
#for ax in axs.flat:
 #   ax.set(xlabel='x-label', ylabel='y-label')

"""


plt.show()
