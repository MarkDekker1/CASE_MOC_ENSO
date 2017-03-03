# ---------------------------------------------- #
# COUPLED SYSTEM: MOC => ENSO
# ---------------------------------------------- #

# ---------------------------------------------- #
# PREAMBULE
# ---------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------- #
# INTEGRATION PARAMETERS
# ---------------------------------------------- #

tmax=400000.
dt=1

# ---------------------------------------------- #
# PHYSICAL PARAMETERS (TAB 10.3, NCD) (IN DAYS)
# ---------------------------------------------- #

#Leading:
tr = 2500.                  # temperature relaxation tscale (NORMALLY 25!)
H = 4500.                   # mean ocean depth
td = 180.*365.              # diffusion time scale
ta = 29.*365.               # advective time scale
q = 192e10*3600.*24.        # transport coefficient
V = 300*4.5*8200.*1e9       # Ocean volume
alphat = 1e-4               # thermal expansion coefficient
alphas = 76e-5              # haline contraction coefficient
S0 = 35                     # Reference salinity
Th = 25                     # meridional temp diff
rho0 = 1.                   # Reference density
F0 = 2.3/365                # Scale freshwater flux
Fs = 0.001
gamma=10

#Following:
Tr0 = 16.
alpha = 1./180.
Hf = 100.
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
        
DT0 = 0.2
DS0 = 0.2
 
T10 = 0
T20 = 0
h10 = 0

# ---------------------------------------------- #
# FORCING AND COUPLING
# ---------------------------------------------- #

def force():
    global Fs
    Fs=min(0.006+max(0.00000005*(t-40000),0),0.0082)
    
def couple(DT):
    global tau_ext
    tau_ext=-0.01+0.01*(DT-23)
    
# ---------------------------------------------- #
# HELPER EQUATIONS
# ---------------------------------------------- #
    
def flow(DT,DS):
    return gamma*Drhof(DT,DS)/rho0

def Q(Drho):
    return 1./td+q*Drho**2./(V)
    
def Drhof(DT,DS):
    return rho0*(alphat*DT-alphas*DS)
    
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

def dDT(t,DT,DS):
    return -1./tr*(DT-Th)-Q(Drhof(DT,DS))*DT

def dDS(t,DT,DS):
    return Fs*S0/H-Q(Drhof(DT,DS))*DS

def dT1(t,T1,T2,h1):
    return -alpha*(T1-Tr)-(T2-T1)*u(T1,T2,h1)/(L/2.)

def dT2(t,T1,T2,h1):
    Tsub = Tr-(Tr-Tr0)/2.*(1.-np.tanh((Hf-z0+h2(T1,T2,h1))/hstar))
    return -alpha*(T2-Tr)-(T2-Tsub)*w(T1,T2,h1)/Hm
    
def dh1(t,T1,T2,h1):
    return r*(-h1-b*L*tau(T1,T2,h1)/2.)

# ---------------------------------------------- #
# TIME INTEGRATION LOOP
# ---------------------------------------------- #

DTvec=[DT0]
DSvec=[DS0]
PSIvec=[DS0]
Fvec=[Fs]
PSIvec=[flow(DT0,DS0)]
T1vec=[T10]
T2vec=[T20]
Hvec=[h10]
tvec=[0]
i=0
while i<tmax/dt:
    t=i*dt  
    force()
    
    k1 = dt*dDT(t,DTvec[i],DSvec[i])
    k2 = dt*dDT(t+0.5*dt,DTvec[i]+k1*0.5,DSvec[i])
    k3 = dt*dDT(t+0.5*dt,DTvec[i]+k2*0.5,DSvec[i])
    k4 = dt*dDT(t+dt,DTvec[i]+k3,DSvec[i])
    DTnew=DTvec[i]+(k1+2.*k2+2.*k3+k4)/6.
    
    k1 = dt*dDS(t,DTvec[i],DSvec[i])
    k2 = dt*dDS(t+0.5*dt,DTvec[i],DSvec[i]+k1*0.5)
    k3 = dt*dDS(t+0.5*dt,DTvec[i],DSvec[i]+k2*0.5)
    k4 = dt*dDS(t+dt,DTvec[i],DSvec[i]+k3)
    DSnew=DSvec[i]+(k1+2.*k2+2.*k3+k4)/6.
    
    DTvec.append(DTnew)
    DSvec.append(DSnew)
    
    Fvec.append(Fs)
    PSIvec.append(flow(DTnew,DSnew))
    
    couple(DTnew)
    
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
    i=i+1

tvec=np.array(tvec)/1.
#%%
fig, ax1 = plt.subplots(figsize=(12,4))
ax1.plot(tvec,np.array(Fvec),'k--',linewidth=3)
ax1.plot(tvec,np.array(PSIvec),'--',linewidth=3,color='Grey')
ax1.set_ylabel(r'Leading variables',fontsize=15)
ax1.tick_params(axis='both',which='major',labelsize=15)
l1=ax1.legend([r'$F_s$',r'$\Psi$'],loc='upper left')
l1.set_zorder(100)

ax2 = ax1.twinx()
ax2.plot(tvec,T1vec,'b',linewidth=3)
ax2.plot(tvec,T2vec,'r',linewidth=3)
ax2.plot(tvec,np.array(Hvec),'y',linewidth=3)
ax2.set_xlim([0,tvec[-1]])
#plt.ylim([-0.002,0.02])
l2=ax2.legend([r'$T_1$',r'$T_2$',r'$h_1$'],loc='upper right')
l2.set_zorder(10)
ax2.set_xlabel('Time',fontsize=15)
ax2.set_ylabel(r'Following variables',fontsize=15)
ax2.tick_params(axis='both',which='major',labelsize=15)