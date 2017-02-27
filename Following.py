# ---------------------------------------------- #
# LEADING SYSTEM: ENSO => ROBERTS ET AL. 2016
# ---------------------------------------------- #

# Info gained from Robert et al. 2016 and NCD, Ch 8, p 46
# Still dimensional version used
# I couple this system with the leading system
# using the function text(DT)

# ---------------------------------------------- #
# PREAMBULE
# ---------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------- #
# INTEGRATION PARAMETERS
# ---------------------------------------------- #

tmax=5000.
dt=0.1

# ---------------------------------------------- #
# PHYSICAL PARAMETERS (TAB 8.6, NCD, TAB 1, Roberts)
# ---------------------------------------------- #

Tr0 = 16
Tr = 29.5
alpha = 1/180.
H = 100
hstar = 62
eps = 0.11
L = 15e6
r = 1/400.
z0 = 75
mu = 0.0026
zeta = 1.3
kwant = 22
        
# ---------------------------------------------- #
# INITIAL CONDITIONS
# ---------------------------------------------- #
        
T10 = 1
T20 = 1
h10 = 1

# ---------------------------------------------- #
# HELPER EQUATIONS
# ---------------------------------------------- #

def Tsub(T1,T2,h1):
    return (Tr+Tr0)/2.-(Tr-Tr0)/2.*np.tanh((H-z0+h1+kwant*(T2-T1))/hstar)
    
def tau(T1,T2):
    return tau_ext(25)+mu*(T2-T1)
    
def tau_ext(DT):
    return DT/100.*mu

# ---------------------------------------------- #
# DIFFERENTIAL EQUATIONS
# ---------------------------------------------- #

def dT1(DT,T1,T2,h1):
    return -alpha*(T1-Tr)-eps*tau(T1,T2)*(T2-T1)

def dT2(DT,T1,T2,h1):
    return -alpha*(T2-Tr)-zeta*tau(T1,T2)*(T2-Tsub(T1,T2,h1))
    
def dh1(DT,T1,T2,h1):
    return r*(-h1-kwant*(T2-T1)/2.)

# ---------------------------------------------- #
# TIME INTEGRATION LOOP
# ---------------------------------------------- #

T1vec=[T10]
T2vec=[T20]
Hvec=[h10]
Cvec=[tau(T10,T20)]
tvec=[0]
i=0
while i<tmax/dt:
    t=i*dt
    
    k1 = dt*dT1(t,T1vec[i],T2vec[i],Hvec[i])
    k2 = dt*dT1(t+0.5*dt,T1vec[i]+k1*0.5,T2vec[i],Hvec[i])
    k3 = dt*dT1(t+0.5*dt,T1vec[i]+k2*0.5,T2vec[i],Hvec[i])
    k4 = dt*dT1(t+dt,T1vec[i]+k3,T2vec[i],Hvec[i])
    T1new=T1vec[i]+(k1+2.*k2+2.*k3+k4)/6.
    
    k1 = dt*dT2(t,T1vec[i],T2vec[i],Hvec[i])
    k2 = dt*dT2(t+0.5*dt,T1vec[i],T2vec[i]+k1*0.5,Hvec[i])
    k3 = dt*dT2(t+0.5*dt,T1vec[i],T2vec[i]+k2*0.5,Hvec[i])
    k4 = dt*dT2(t+dt,T1vec[i],T2vec[i]+k3,Hvec[i])
    T2new=T2vec[i]+(k1+2.*k2+2.*k3+k4)/6.
    
    k1 = dt*dh1(t,T1vec[i],T2vec[i],Hvec[i])
    k2 = dt*dh1(t+0.5*dt,T1vec[i],T2vec[i],Hvec[i]+k1*0.5)
    k3 = dt*dh1(t+0.5*dt,T1vec[i],T2vec[i],Hvec[i]+k2*0.5)
    k4 = dt*dh1(t+dt,T1vec[i],T2vec[i],Hvec[i]+k3)
    h1new=Hvec[i]+(k1+2.*k2+2.*k3+k4)/6.
    
    T1vec.append(T1new)
    T2vec.append(T2new)
    Hvec.append(h1new)
    
    Cvec.append(tau(T1new,T2new))
    tvec.append(t)
    i=i+1

tvec=np.array(tvec)/365.

plt.plot(tvec,T1vec,linewidth=3)
plt.plot(tvec,T2vec,linewidth=3)
plt.plot(tvec,np.array(Hvec),linewidth=3)
plt.legend([r'$T_1$',r'$T_2$',r'$h_1$'],loc='best')
plt.xlim([0,tvec[-1]])
#plt.ylim([0,30])