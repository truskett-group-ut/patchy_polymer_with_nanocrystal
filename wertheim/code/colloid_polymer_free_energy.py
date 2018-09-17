from __future__ import division
import numpy as np

#free energy per V from assembling the polymer mixture
def Acp(p, x_c, d_c, d_p, m):
    
    #just convert the notation to what was used in mathematica
    xc = x_c
    dc = d_c
    dp = d_p
    
    #yup, this is horrible
    a = (24*np.power(dp,9)*np.power(m,4)*p*np.pi - 3*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2) + 
     42*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*xc + 36*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*xc + 
     18*dc*np.power(dp,8)*np.power(m,3)*p*np.pi*xc - 96*np.power(dp,9)*np.power(m,4)*p*np.pi*xc - 
     9*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc - 
     3*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc - 
     3*dc*np.power(dp,11)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc + 
     15*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*xc + 
     18*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,2) + 
     54*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,2) + 
     54*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,2) + 
     18*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,2) - 
     126*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,2) - 
     108*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,2) - 54*dc*np.power(dp,8)*np.power(m,3)*p*np.pi*np.power(xc,2) + 
     144*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,2) - 
     9*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) - 
     9*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) - 
     9*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) - 
     3*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) + 
     36*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) + 
     12*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) + 
     12*dc*np.power(dp,11)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) - 
     30*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2) + 18*np.power(dc,8)*dp*m*p*np.pi*np.power(xc,3) + 
     36*np.power(dc,7)*np.power(dp,2)*m*p*np.pi*np.power(xc,3) + 42*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,3) - 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,3) - 
     108*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,3) - 
     108*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,3) - 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,3) + 
     126*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,3) + 
     108*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,3) + 54*dc*np.power(dp,8)*np.power(m,3)*p*np.pi*np.power(xc,3) - 
     96*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,3) - 
     3*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) - 
     9*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) - 
     9*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) - 
     9*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) + 
     27*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) + 
     27*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) + 
     27*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) + 
     9*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) - 
     54*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) - 
     18*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) - 
     18*dc*np.power(dp,11)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) + 
     30*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3) + 24*np.power(dc,9)*p*np.pi*np.power(xc,4) - 
     18*np.power(dc,8)*dp*m*p*np.pi*np.power(xc,4) - 36*np.power(dc,7)*np.power(dp,2)*m*p*np.pi*np.power(xc,4) - 
     42*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,4) + 18*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,4) + 
     54*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,4) + 
     54*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,4) + 
     18*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,4) - 
     42*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,4) - 
     36*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,4) - 18*dc*np.power(dp,8)*np.power(m,3)*p*np.pi*np.power(xc,4) + 
     24*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,4) - 3*np.power(dc,11)*dp*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     3*np.power(dc,10)*np.power(dp,2)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     9*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) + 
     6*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) + 
     18*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) + 
     18*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) + 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     27*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     27*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     27*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     9*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) + 
     36*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) + 
     12*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) + 
     12*dc*np.power(dp,11)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     15*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4) - 
     3*np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 3*np.power(dc,11)*dp*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 
     3*np.power(dc,10)*np.power(dp,2)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 
     9*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) - 
     3*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) - 
     9*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) - 
     9*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) - 
     9*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 
     9*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 
     9*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 
     9*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 
     3*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) - 
     9*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) - 
     3*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) - 
     3*dc*np.power(dp,11)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 
     3*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5) + 72*np.power(dp,6)*np.power(m,2)*np.log(6) - 
     72*np.power(dp,6)*np.power(m,3)*np.log(6) - 24*np.power(dp,9)*np.power(m,3)*p*np.pi*np.log(6) + 
     24*np.power(dp,9)*np.power(m,4)*p*np.pi*np.log(6) + 2*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.log(6) - 
     2*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.log(6) + 144*np.power(dc,3)*np.power(dp,3)*m*xc*np.log(6) - 
     72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*xc*np.log(6) - 108*np.power(dc,2)*np.power(dp,4)*np.power(m,2)*xc*np.log(6) - 
     180*np.power(dp,6)*np.power(m,2)*xc*np.log(6) + 216*np.power(dp,6)*np.power(m,3)*xc*np.log(6) - 
     72*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*xc*np.log(6) + 48*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*xc*np.log(6) + 
     36*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*xc*np.log(6) + 84*np.power(dp,9)*np.power(m,3)*p*np.pi*xc*np.log(6) - 
     96*np.power(dp,9)*np.power(m,4)*p*np.pi*xc*np.log(6) + 
     8*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(6) - 
     6*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(6) - 
     3*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(6) - 
     9*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(6) + 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(6) + 72*np.power(dc,6)*np.power(xc,2)*np.log(6) - 
     36*np.power(dc,6)*m*np.power(xc,2)*np.log(6) - 108*np.power(dc,4)*np.power(dp,2)*m*np.power(xc,2)*np.log(6) - 
     216*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,2)*np.log(6) + 144*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,2)*np.log(6) + 
     216*np.power(dc,2)*np.power(dp,4)*np.power(m,2)*np.power(xc,2)*np.log(6) + 144*np.power(dp,6)*np.power(m,2)*np.power(xc,2)*np.log(6) - 
     216*np.power(dp,6)*np.power(m,3)*np.power(xc,2)*np.log(6) - 72*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,2)*np.log(6) + 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,2)*np.log(6) + 
     36*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,2)*np.log(6) + 
     36*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,2)*np.log(6) + 
     180*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,2)*np.log(6) - 
     144*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,2)*np.log(6) - 
     108*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,2)*np.log(6) - 
     108*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,2)*np.log(6) + 144*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,2)*np.log(6) + 
     12*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) - 
     7*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) - 
     6*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) - 
     3*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) - 
     28*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) + 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) + 
     12*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) + 
     16*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) - 
     20*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(6) - 72*np.power(dc,6)*np.power(xc,3)*np.log(6) + 
     36*np.power(dc,6)*m*np.power(xc,3)*np.log(6) + 108*np.power(dc,4)*np.power(dp,2)*m*np.power(xc,3)*np.log(6) + 
     72*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,3)*np.log(6) - 72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,3)*np.log(6) - 
     108*np.power(dc,2)*np.power(dp,4)*np.power(m,2)*np.power(xc,3)*np.log(6) - 36*np.power(dp,6)*np.power(m,2)*np.power(xc,3)*np.log(6) + 
     72*np.power(dp,6)*np.power(m,3)*np.power(xc,3)*np.log(6) - 24*np.power(dc,9)*p*np.pi*np.power(xc,3)*np.log(6) + 
     12*np.power(dc,9)*m*p*np.pi*np.power(xc,3)*np.log(6) + 36*np.power(dc,7)*np.power(dp,2)*m*p*np.pi*np.power(xc,3)*np.log(6) + 
     120*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,3)*np.log(6) - 
     72*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,3)*np.log(6) - 
     72*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,3)*np.log(6) - 
     72*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,3)*np.log(6) - 
     144*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,3)*np.log(6) + 
     144*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,3)*np.log(6) + 
     108*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,3)*np.log(6) + 
     60*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,3)*np.log(6) - 96*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,3)*np.log(6) + 
     8*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) - 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) - 
     3*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) - 
     6*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) - 
     31*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) + 
     21*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) + 
     18*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) + 
     9*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) + 
     36*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) - 
     36*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) - 
     18*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) - 
     14*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) + 
     20*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(6) + 
     24*np.power(dc,9)*p*np.pi*np.power(xc,4)*np.log(6) - 12*np.power(dc,9)*m*p*np.pi*np.power(xc,4)*np.log(6) - 
     36*np.power(dc,7)*np.power(dp,2)*m*p*np.pi*np.power(xc,4)*np.log(6) - 48*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,4)*np.log(6) + 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(6) + 
     36*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(6) + 
     36*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(6) + 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(6) - 
     48*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,4)*np.log(6) - 
     36*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,4)*np.log(6) - 
     12*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,4)*np.log(6) + 24*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,4)*np.log(6) + 
     2*np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     3*np.power(dc,10)*np.power(dp,2)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     14*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) + 
     8*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) + 
     6*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) + 
     12*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) + 
     26*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     21*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     18*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     9*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     20*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) + 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) + 
     12*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) + 
     6*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6) - 
     2*np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     3*np.power(dc,10)*np.power(dp,2)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     6*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 
     3*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 
     6*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 
     7*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     7*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     6*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     3*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 
     6*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 
     3*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 
     np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) + 
     2*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6) - 36*np.power(dp,6)*np.power(m,2)*np.log(12) + 
     36*np.power(dp,6)*np.power(m,3)*np.log(12) + 12*np.power(dp,9)*np.power(m,3)*p*np.pi*np.log(12) - 
     12*np.power(dp,9)*np.power(m,4)*p*np.pi*np.log(12) - np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.log(12) + 
     np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.log(12) - 72*np.power(dc,3)*np.power(dp,3)*m*xc*np.log(12) + 
     72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*xc*np.log(12) + 108*np.power(dp,6)*np.power(m,2)*xc*np.log(12) - 
     108*np.power(dp,6)*np.power(m,3)*xc*np.log(12) + 36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*xc*np.log(12) - 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*xc*np.log(12) - 48*np.power(dp,9)*np.power(m,3)*p*np.pi*xc*np.log(12) + 
     48*np.power(dp,9)*np.power(m,4)*p*np.pi*xc*np.log(12) - 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(12) + 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(12) + 
     5*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(12) - 
     5*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(12) - 36*np.power(dc,6)*np.power(xc,2)*np.log(12) + 
     36*np.power(dc,6)*m*np.power(xc,2)*np.log(12) + 144*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,2)*np.log(12) - 
     144*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,2)*np.log(12) - 108*np.power(dp,6)*np.power(m,2)*np.power(xc,2)*np.log(12) + 
     108*np.power(dp,6)*np.power(m,3)*np.power(xc,2)*np.log(12) + 36*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,2)*np.log(12) - 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,2)*np.log(12) - 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,2)*np.log(12) + 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,2)*np.log(12) + 
     72*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,2)*np.log(12) - 72*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,2)*np.log(12) - 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(12) + 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(12) + 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(12) - 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(12) - 
     10*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(12) + 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*np.log(12) + 36*np.power(dc,6)*np.power(xc,3)*np.log(12) - 
     36*np.power(dc,6)*m*np.power(xc,3)*np.log(12) - 72*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,3)*np.log(12) + 
     72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,3)*np.log(12) + 36*np.power(dp,6)*np.power(m,2)*np.power(xc,3)*np.log(12) - 
     36*np.power(dp,6)*np.power(m,3)*np.power(xc,3)*np.log(12) + 12*np.power(dc,9)*p*np.pi*np.power(xc,3)*np.log(12) - 
     12*np.power(dc,9)*m*p*np.pi*np.power(xc,3)*np.log(12) - 72*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,3)*np.log(12) + 
     72*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,3)*np.log(12) + 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,3)*np.log(12) - 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,3)*np.log(12) - 
     48*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,3)*np.log(12) + 48*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,3)*np.log(12) - 
     4*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) + 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) + 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) - 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) - 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) + 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) + 
     10*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) - 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*np.log(12) - 
     12*np.power(dc,9)*p*np.pi*np.power(xc,4)*np.log(12) + 12*np.power(dc,9)*m*p*np.pi*np.power(xc,4)*np.log(12) + 
     36*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,4)*np.log(12) - 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(12) - 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(12) + 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,4)*np.log(12) + 
     12*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,4)*np.log(12) - 12*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,4)*np.log(12) - 
     np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) + 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) + 
     8*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) - 
     8*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) - 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) + 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) + 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) - 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) - 
     5*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) + 
     5*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(12) + 
     np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) - 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) - 
     4*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) + 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) + 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) - 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) - 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) + 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) + 
     np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) - 
     np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(12) + 
     36*np.power(dp,6)*np.power(m,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dp,6)*np.power(m,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     12*np.power(dp,9)*np.power(m,3)*p*np.pi*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     12*np.power(dp,9)*np.power(m,4)*p*np.pi*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     72*np.power(dc,3)*np.power(dp,3)*m*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     108*np.power(dp,6)*np.power(m,2)*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     108*np.power(dp,6)*np.power(m,3)*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     48*np.power(dp,9)*np.power(m,3)*p*np.pi*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     48*np.power(dp,9)*np.power(m,4)*p*np.pi*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     5*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     5*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*xc*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     36*np.power(dc,6)*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dc,6)*m*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     144*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     144*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     108*np.power(dp,6)*np.power(m,2)*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     108*np.power(dp,6)*np.power(m,3)*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     72*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     72*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,2)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     10*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dc,6)*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     36*np.power(dc,6)*m*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     72*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dp,6)*np.power(m,2)*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     36*np.power(dp,6)*np.power(m,3)*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     12*np.power(dc,9)*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     12*np.power(dc,9)*m*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     72*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     72*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     48*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     48*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     4*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     10*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     12*np.power(dc,9)*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     12*np.power(dc,9)*m*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     12*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     12*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     8*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     8*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     5*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     5*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     4*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) + 
     np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc) - 
     108*np.power(dp,6)*np.power(m,2)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dp,6)*np.power(m,3)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dp,9)*np.power(m,3)*p*np.pi*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dp,9)*np.power(m,4)*p*np.pi*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     3*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     216*np.power(dc,3)*np.power(dp,3)*m*xc*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     144*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,2)*np.power(dp,4)*np.power(m,2)*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     288*np.power(dp,6)*np.power(m,2)*xc*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     324*np.power(dp,6)*np.power(m,3)*xc*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     84*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     132*np.power(dp,9)*np.power(m,3)*p*np.pi*xc*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     144*np.power(dp,9)*np.power(m,4)*p*np.pi*xc*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     10*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     14*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     15*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dc,6)*np.power(xc,2)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dc,6)*m*np.power(xc,2)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,4)*np.power(dp,2)*m*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     360*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     288*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     216*np.power(dc,2)*np.power(dp,4)*np.power(m,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     252*np.power(dp,6)*np.power(m,2)*np.power(xc,2)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) + 324*np.power(dp,6)*np.power(m,3)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     72*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     288*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     252*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     180*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     216*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     13*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     6*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     44*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     40*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     26*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     30*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,6)*np.power(xc,3)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     72*np.power(dc,6)*m*np.power(xc,3)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dc,4)*np.power(dp,2)*m*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     144*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     144*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,2)*np.power(dp,4)*np.power(m,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dp,6)*np.power(m,2)*np.power(xc,3)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) - 108*np.power(dp,6)*np.power(m,3)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,9)*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     24*np.power(dc,9)*m*p*np.pi*np.power(xc,3)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,7)*np.power(dp,2)*m*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     192*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     144*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     252*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     252*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     144*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     8*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     6*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     49*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     39*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     18*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     9*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     60*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     60*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     18*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     24*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     30*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,9)*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     24*np.power(dc,9)*m*p*np.pi*np.power(xc,4)*np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,7)*np.power(dp,2)*m*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     84*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     72*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,5)*np.power(dp,4)*np.power(m,2)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,4)*np.power(dp,5)*np.power(m,2)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     72*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     84*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,2)*np.power(dp,7)*np.power(m,3)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     24*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     3*np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     2*np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dc,10)*np.power(dp,2)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     22*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     16*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     6*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     44*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     39*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     18*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     9*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     40*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     11*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     15*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     2*np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     3*np.power(dc,10)*np.power(dp,2)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     10*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     8*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dc,8)*np.power(dp,4)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     6*np.power(dc,7)*np.power(dp,5)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     13*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     13*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     6*np.power(dc,5)*np.power(dp,7)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     3*np.power(dc,4)*np.power(dp,8)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     8*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     10*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     3*np.power(dc,2)*np.power(dp,10)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     2*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     3*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dp,6)*np.power(m,2)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) - 36*np.power(dp,6)*np.power(m,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dp,9)*np.power(m,3)*p*np.pi*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) + 12*np.power(dp,9)*np.power(m,4)*p*np.pi*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dc,3)*np.power(dp,3)*m*xc*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) - 72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dp,6)*np.power(m,2)*xc*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) + 108*np.power(dp,6)*np.power(m,3)*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     48*np.power(dp,9)*np.power(m,3)*p*np.pi*xc*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + 
        np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     48*np.power(dp,9)*np.power(m,4)*p*np.pi*xc*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + 
        np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     5*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     5*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*xc*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,6)*np.power(xc,2)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) - 36*np.power(dc,6)*m*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     144*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     144*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dp,6)*np.power(m,2)*np.power(xc,2)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + 
        np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dp,6)*np.power(m,3)*np.power(xc,2)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + 
        np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     72*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     10*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,2)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,6)*np.power(xc,3)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) + 36*np.power(dc,6)*m*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dc,3)*np.power(dp,3)*m*np.power(xc,3)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + 
        np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     72*np.power(dc,3)*np.power(dp,3)*np.power(m,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dp,6)*np.power(m,2)*np.power(xc,3)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + 
        np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dp,6)*np.power(m,3)*np.power(xc,3)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + 
        np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dc,9)*p*np.pi*np.power(xc,3)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) + 12*np.power(dc,9)*m*p*np.pi*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     72*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     72*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     108*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     48*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     48*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     4*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     24*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     10*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     10*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,3)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     12*np.power(dc,9)*p*np.pi*np.power(xc,4)*np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + 
        np.power(dp,3)*m*p*np.pi*xc) - 12*np.power(dc,9)*m*p*np.pi*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,6)*np.power(dp,3)*m*p*np.pi*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,6)*np.power(dp,3)*np.power(m,2)*p*np.pi*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,2)*p*np.pi*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     36*np.power(dc,3)*np.power(dp,6)*np.power(m,3)*p*np.pi*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     12*np.power(dp,9)*np.power(m,3)*p*np.pi*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     12*np.power(dp,9)*np.power(m,4)*p*np.pi*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     8*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     8*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     18*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     16*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     5*np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     5*np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,4)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     np.power(dc,12)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     np.power(dc,12)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     4*np.power(dc,9)*np.power(dp,3)*m*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     4*np.power(dc,9)*np.power(dp,3)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,2)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     6*np.power(dc,6)*np.power(dp,6)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,3)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     4*np.power(dc,3)*np.power(dp,9)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) - 
     np.power(dp,12)*np.power(m,4)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc) + 
     np.power(dp,12)*np.power(m,5)*np.power(p,2)*np.power(np.pi,2)*np.power(xc,5)*
      np.log(12 - np.power(dp,3)*m*p*np.pi - 2*np.power(dc,3)*p*np.pi*xc + np.power(dc,2)*dp*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc))/(np.power(-(np.power(dp,3)*m) - np.power(dc,3)*xc + np.power(dp,3)*m*xc,2)*
     np.power(6 - np.power(dp,3)*m*p*np.pi - np.power(dc,3)*p*np.pi*xc + np.power(dp,3)*m*p*np.pi*xc,2))
    
    #return the free energy per volume, 
    #not per number of sites, 
    #as the formula calculates
    return a*p