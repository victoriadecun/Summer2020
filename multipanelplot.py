import matplotlib.pyplot as plt
import numpy as np

#we skip ruth 9.4 because LISA won't observe it.

x,y = np.loadtxt('ruth_1.75_data.data', delimiter=',', unpack=True)
xfilt,yfilt = np.loadtxt('ruth_1.75_data_filt.data', delimiter=',', unpack=True)

x2,y2 = np.loadtxt('ruth_3.7_data.data', delimiter=',', unpack=True)
x2filt,y2filt = np.loadtxt('ruth_3.7_data_filt.data', delimiter=',', unpack=True)

x3,y3 = np.loadtxt('sandra_0.39_data.data', delimiter=',', unpack=True)
x3filt,y3filt = np.loadtxt('sandra_0.39_data_filt.data', delimiter=',', unpack=True)

x4,y4 = np.loadtxt('sonia_4.6_data.data', delimiter=',', unpack=True)
x4filt,y4filt = np.loadtxt('sonia_4.6_data_filt.data', delimiter=',', unpack=True)


fig, axs = plt.subplots(2,2,figsize=(8,7))

ax=fig.add_subplot(111,frameon=False)
ax.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
ax.set_xlabel('log R [kpc]', fontsize=14, labelpad=15)
ax.set_ylabel(r'log $\rho_{\star}$ [M$_{\odot}$ kpc$^{-3}$]', fontsize=14, labelpad=15)

#'tab:blue' etc deleted

axs[0, 0].plot(x, y, color='gray')
axs[0,0].text(-1.8,6.3,  'z = 1.75', fontsize=12, bbox={'facecolor':'white', 'alpha': 0.5, 'pad': 10})
axs[0,0].plot(xfilt,yfilt,'.', color='blue')
m,b = np.polyfit(xfilt,yfilt,1)
axs[0,0].plot(xfilt,m*xfilt+b,color='red')


axs[0, 1].plot(x2, y2,  color='gray')
axs[0,1].text(-1.8,5.2, 'z = 3.79', fontsize=12, bbox={'facecolor':'white', 'alpha': 0.5, 'pad': 10})
axs[0,1].plot(x2filt,y2filt,'.', color='blue')
m,b = np.polyfit(x2filt,y2filt,1)
axs[0,1].plot(x2filt,m*x2filt+b,color='red')


axs[1, 0].plot(x3, y3, color='gray')
axs[1,0].text(-1.6,6.9, 'z = 0.392', fontsize=12, bbox={'facecolor':'white', 'alpha': 0.5, 'pad': 10})
axs[1,0].plot(x3filt,y3filt,'.', color='blue')
m,b = np.polyfit(x3filt,y3filt,1)
axs[1,0].plot(x3filt,m*x3filt+b,color='red')


axs[1, 1].plot(x4, y4,color='gray')
axs[1,1].text(-1.7,5.1, 'z = 4.60', fontsize=12,bbox={'facecolor':'white', 'alpha': 0.5, '\
pad': 10})
axs[1,1].plot(x4filt,y4filt,'.', color='blue')
m,b = np.polyfit(x4filt,y4filt,1)
axs[1,1].plot(x4filt,m*x4filt+b,color='red')



plt.show()
