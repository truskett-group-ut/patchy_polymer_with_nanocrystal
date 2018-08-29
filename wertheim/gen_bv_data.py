import os, sys
import numpy as np
from scipy import integrate

#calculate the average mayer function for werthim theory
#between off-centered the patch on a particle and the centered
#patch on another particle using a lammps-style gaussian interaction
def Favg(r, rdis, A, B, rcut):  
    
    #place polymer end along the z axis
    polymer_end = np.array([0, 0, r])
    
    #interaction potential
    def SinFpp(theta, phi):
        
        patch = np.array([np.sin(theta)*np.cos(phi),
                          np.sin(theta)*np.sin(phi), 
                          np.cos(theta)])
        patch = patch*rdis
        rpp = np.linalg.norm(polymer_end - patch)
        upp = -A*np.exp(-B*rpp**2)*np.heaviside(rcut - rpp, 0)
        fpp = np.exp(-upp) - 1.0
        
        return np.sin(theta)*fpp
    
    #perform the 2D integration
    integral = integrate.dblquad(SinFpp, 0, 2.0*np.pi, lambda phi: 0, lambda phi: np.pi)
    average = integral[0] / (4.0*np.pi)
    
    return average
    
dp = 1.0
dc = 5.0

Amin = 0.0
Amax = 20.0
Ainc = 0.02
B = 1.0/0.2**2

rdis = dp*2.0**(1.0/6.0) + (dp + dc)/2.0 - dp
rcut = 0.5*dp
rmin = rdis - rcut*1.1
rmax = rdis + rcut*1.1

dr = 0.005
rs = np.arange(rmin, rmax, dr)
As = np.arange(Amin, Amax, Ainc)

path = './_favg__dp={}__dc={}__B={}__rdis={}__rcut={}'.format(dp, dc, B, rdis, rcut)
if not os.path.exists(path):
    os.makedirs(path)
    
for A in As:
    print "Working on A={}".format(A)
    
    data = [A]
    for r in rs:
        data.append(Favg(r, rdis, A, B, rcut))
    
    #convert to an array and save
    data = np.array(data)
    np.savetxt('{}/A={}.txt'.format(path, A), data)