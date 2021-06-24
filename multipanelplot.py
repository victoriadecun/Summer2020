import matplotlib.pyplot as plt
import numpy as np

#we skip ruth 9.4 because LISA won't observe it.
x,y = np.loadtxt('ruth_1.75_data.data', delimiter=',', unpack=True)
x2,y2 = np.loadtxt('ruth_3.7_data.data', delimiter=',', unpack=True)
x3,y3 = np.loadtxt('sandra_0.39_data.data', delimiter=',', unpack=True)
x4,y4 = np.loadtxt('sonia_4.6_data.data', delimiter=',', unpack=True)


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
axs[0,1].text(-1.45,7.60, 'z = 3.79', fontsize=12)
axs[0,1].set_xlim([-1.5,-0.3])
axs[0,1].set_ylim([7.50,8.75])
m,b = np.polyfit(x2,y2,1)
axs[0,1].plot(x2,m*x2+b,color='black')


axs[1, 0].plot(x3, y3, 'tab:green')
#axs[1, 0].set_title('Axis [1, 0]')
axs[1,0].set_ylabel(r'log $\rho_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
axs[1,0].set_xlabel('log R [kpc]')
axs[1,0].text(-1.4,7.6, 'z = 0.392', fontsize=12)
axs[1,0].set_xlim([-1.5,-0.3])
axs[1,0].set_ylim([7.50,8.25])
m,b = np.polyfit(x3,y3,1)
axs[1,0].plot(x3,m*x3+b,color='black')


axs[1, 1].plot(x4, y4, 'tab:red')
#axs[1, 1].set_title('Axis [1, 1]')
axs[1,1].set_xlabel('log R [kpc]')
axs[1,1].text(-1.4,7.8, 'z = 4.60', fontsize=12)
m,b = np.polyfit(x4,y4,1)
axs[1,1].plot(x4,m*x4+b,color='black')


#for ax in axs.flat:
 #   ax.set(xlabel='x-label', ylabel='y-label')




plt.show()
