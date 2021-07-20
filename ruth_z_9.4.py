import pynbody
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from matplotlib.ticker import NullFormatter

#load snapshot

s=pynbody.load("/mnt/data0/jillian/h229/h229.cosmo50PLK.3072gst5HbwK1BH.000139/h229.cosmo50PLK.3072gst5HbwK1BH.000139")

#convert the units
s.physical_units()

"""
#  load any available halo
h = s.halos()
# function to find black hole
def findBH(s):
    BHfilter = pynbody.filt.LowPass('tform',0.0)
    BH = s.stars[BHfilter]
    return BH
BH = findBH(s)
print BH
#function to find the halos that the galaxy is in
def findBHhalos(s):
    BH = findBH(s)
    BHhalos = BH['amiga.grp']
    return BHhalos
#using the function the halos
BHhalos = findBHhalos(s)
#printing the halos
print BHhalos
#sorting the halos, indexes/indecis are like an exact address
currenthalo = np.argsort(BHhalos)
print BHhalos[currenthalo]
def getz(s):
    return s.properties['z']
def gettime(s):
    return pynbody.analysis.cosmology.age(s)
for i in currenthalo:
    #which halo are we on?
    currenthalo = BHhalos[i]
    print 'current halo: ', currenthalo
    print i
    
    #put the galaxy you care about in the center of the simulation
    pynbody.analysis.angmom.faceon(h[currenthalo])
    #with pynbody.analysis.halo.center(h[currenthalo], mode='hyb'):
#for halo in BHhalos:
#    halo/galaxy we need to look at is no. 2
#    pynbody.analysis.angmom.faceon(h[halo])
#    prints the black hole
#    print BH
    #this is the position of black hole
    BHposition=BH['pos']
    #printing the position of black hole
    #print BHposition
    #putting the x-values into a column
    BHx= BHposition[[i],0]
    print "x postion", BHx
    #putting the y-values into a column
    BHy= BHposition[[i],1]
    print "y position", BHy
    #putting the z-values into a column
    BHz= BHposition[[i],2]
    print "z position", BHz
    #the .5 is the square root , this is the distance formula
    #distance =((BHx**2)+(BHy**2)+(BHz**2))**(.5)
    #print 'this is the distance :'
    #print "this is the distance :", distance
    data = [currenthalo, BH['iord'][i]] 
    print "This is current halo, ID", data
"""

#  load any available halo                                                                             
h = s.halos()
h4=h[4]

#put galaxy in center of simulation
pynbody.analysis.angmom.faceon(h4)

#create a profile object for stars in 3D                                                               
p = pynbody.analysis.profile.Profile(h[4].s,min=.01,max=2,nbins=50,ndim=3,type='log')

#new r in and r out. calculate virial radius, and 2% of v.r. will be r out.
#original range was out is 0.4, range in 0.04

vr = pynbody.analysis.halo.virial_radius(h4)
z = 9.4247567
rin = (0.68)/(1+z)
rout = 0.02*vr


#filter is filt, range of rbins >in and <out
filt=np.where((p['rbins']>rin) & (p['rbins']<rout))
print(filt)

#show slope of density profile
plt.plot(np.log10(p['rbins']),np.log10(p['density']),'k')
x=np.array(np.log10(p['rbins'][filt]))
y=np.array(np.log10(p['density'][filt]))
m,b=np.polyfit(x,y,1)
plt.plot(x,y,'o')
#plt.xscale('log')
#plt.yscale('log')
#plt.plot.semilogy()
#plt.plot.semilogx()
plt.xlabel('log R [kpc]')
plt.ylabel(r'log $\rho_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
plt.plot(x,m*x+b)
#plt.twinx()
plt.title('Density Profile of Stars')
plt.show()
print ("This is Slope & Y-Intercept: ", m,b)

#save x and y into a data file

#c = np.savetxt('ruth_1.75_data.data', (x,y))
c = np.savetxt('ruth_9.4_data.data', np.column_stack((x,y)), delimiter=',')


#calculate distance for separation
#finds number of BH
def findBH(s):
     BHfilter = np.where((s.stars['iord']==60352791)|(s.stars['iord']==60354626))
     BH = s.stars[BHfilter]
     return BH
BH = findBH(s)
print('number of black holes: ',len(BH))
#how far away BH from galaxy center
#it's x,y,z positions
with pynbody.analysis.halo.center(h[4], mode='hyb'):
    print (h[4]['pos'][0])
    print (h[4]['pos'][1])
    print (h[4]['pos'][2])
print('pos')
#position of BH
BHposition=BH['pos']
print('this is BH position: ', BHposition)
#put x-values in column
BHx1=BHposition[[0],0]
BHx2=BHposition[[1],0]
#print('this is the x-position: ', BHx)
#put y-values in column
BHy1=BHposition[[0],1]
BHy2=BHposition[[1],1]
#print('this is the y-position: ', BHy)
#put z-values in column
BHz1=BHposition[[0],2]
BHz2=BHposition[[1],2]
#print('this is the z-position: ', BHz)
#distance formula between two points
distance=(((BHx2-BHx1)**2)+((BHy2-BHy1)**2)+((BHz2-BHz1)**2))**(0.5)
print('This is the distance between BH: ', distance)

#ignore up to print distance 2
#.5 is the square root; distance
#distance1=((BHx1**2)+(BHy1**2)+(BHz1**2))**(.5)
#distance2=((BHx2**2)+(BHy2**2)+(BHz2**2))**(.5)
#print('this is the distance 1: ', distance1)
#print('this is the distance 2: ', distance2)


#Calculate velocity dispersions
BH_pos=BH['pos']
BHx=BH_pos[:,0]
BHy=BH_pos[:,1]
BHz=BH_pos[:,2]
BH_position=np.array([BHx[0],BHy[0],BHz[0]])
#this puts sphere around BH - 0.116 is how i got 99 stars 
radius= 1.0  #kpc
sphere = pynbody.filt.Sphere(radius, cen =(BH_position))
stars = s.stars[0:]
in_sphere = stars[sphere]
total_stars = len(in_sphere)
print("Total stars: ",total_stars)
#find velocity
velocity=in_sphere['vel']
#find velocity of stars in x,y,z
x = np.std(velocity[0])
y = np.std(velocity[1])
z = np.std(velocity[2])
#average of these by dividing by total (velocity dispersion)
vel_answer = np.sqrt((x)**2 + (y)**2 + (z)**2)
print("Velocity Dispersion: ",vel_answer)


