#!python
from __future__ import division,print_function
from pyPRISM.potential.Potential import Potential
import numpy as np

class ShiftedWeeksChandlerAnderson(Potential):    
 
    def __init__(self,epsilon,length_scale,sigma,r_tolerance=0.01):
        
        self.epsilon = epsilon
        self.length_scale = length_scale
        self.sigma = sigma
        self.delta = self.sigma - self.length_scale
        self.r_tolerance = r_tolerance
        self.ur  = lambda r: 4.0 * self.epsilon * ((self.length_scale/(r-self.delta))**(12.0) - (self.length_scale/(r-self.delta))**(6.0) + 1.0/4.0)
        
    def __repr__(self):
        
        return '<Potential: ShiftedWeeksChandlerAnderson>'
        
    def calculate(self,r):
        
        ur = 0.0*r
        mask = (r-self.delta) > self.r_tolerance 
        ur[mask] = self.ur(r[mask])
        
        mask = (r-self.delta) <= self.r_tolerance
        ur[mask] = ur[np.argmax(ur)]
        
        mask = r>(2.0**(1.0/6.0))*self.length_scale + self.delta
        ur[mask] = 0.0
                
        return ur       
