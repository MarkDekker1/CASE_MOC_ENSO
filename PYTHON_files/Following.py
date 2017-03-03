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

tmax=30000
dt=1

# ---------------------------------------------- #
# PHYSICAL PARAMETERS (TAB 8.6, NCD, TAB 1, Roberts)
# ---------------------------------------------- #

Tr0 = 16.
alpha = 1./180.
H = 100.
hstar = 62.
L = 15e6
r = 1./400.
z0 = 75.
mu = 0.0026
zeta = 1.3
Q = 22.
beta=mu/0.02
b = 22*beta/(mu*L)
Tr = 29.5
Hm = 50.
tau_ext=-0.01
eps=0.1

# ---------------------------------------------- #
# INITIAL CONDITIONS
# ---------------------------------------------- #
        
T10 = 0
T20 = 0
h10 = 0

# ---------------------------------------------- #
# FORCING
# ---------------------------------------------- #

def force():
    global tau_ext
    tau_ext=min(-0.01+max(t-10000,0)*0.00001,0.00001)

# ---------------------------------------------- #
# HELPER EQUATIONS
# ---------------------------------------------- #

def tau(T1,T2,h1):
    return tau_ext+mu/beta*(T2-T1)
    
def u(T1,T2,h1):
    return eps*beta*tau(T1,T2,h1)*L/2.
    
def w(T1,T2,h1):
    return -zeta*beta*tau(T1,T2,h1)*Hm
    
def h2(T1,T2,h1):
    return h1+b*L*tau(T1,T2,h1)

# ---------------------------------------------- #
# DIFFERENTIAL EQUATIONS
# ---------------------------------------------- #

def dT1(t,T1,T2,h1):
    return -alpha*(T1-Tr)-(T2-T1)*u(T1,T2,h1)/(L/2.)

def dT2(t,T1,T2,h1):
    Tsub = Tr-(Tr-Tr0)/2.*(1.-np.tanh((H-z0+h2(T1,T2,h1))/hstar))
    return -alpha*(T2-Tr)-(T2-Tsub)*w(T1,T2,h1)/Hm
    
def dh1(t,T1,T2,h1):
    return r*(-h1-b*L*tau(T1,T2,h1)/2.)

# ---------------------------------------------- #
# TIME INTEGRATION LOOP
# ---------------------------------------------- #

T1vec=[T10]
T2vec=[T20]
Hvec=[h10]
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
    
    tvec.append(t)
    force()
    i=i+1

tvec=np.array(tvec)/1.

plt.plot(tvec,T1vec,linewidth=3)
plt.plot(tvec,T2vec,linewidth=3)
plt.plot(tvec,np.array(Hvec),linewidth=3)
plt.legend([r'$T_1$',r'$T_2$',r'$h_1$'],loc='best')
plt.xlim([5000,tvec[-1]])
plt.xlabel('Time',fontsize=15)
plt.ylabel(r'Variables',fontsize=15)
plt.tick_params(axis='both',which='major',labelsize=15)