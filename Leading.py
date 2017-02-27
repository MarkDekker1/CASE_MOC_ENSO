# ---------------------------------------------- #
# LEADING SYSTEM: MOC => STOMMEL (CH 10, p 29)
# ---------------------------------------------- #

# For now I use the dimensional version of Stommel.
# See page 30 of CH 10 of Nonlinear Climate Dynamics
# for also a nondimensional version.
# I force this system using Fs(t)

# ---------------------------------------------- #
# PREAMBULE
# ---------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------- #
# INTEGRATION PARAMETERS
# ---------------------------------------------- #

tmax=100000.
dt=1

# ---------------------------------------------- #
# PHYSICAL PARAMETERS (TAB 10.3, NCD) (IN DAYS)
# ---------------------------------------------- #

tr = 25.                    # temperature relaxation tscale
H = 4500.                   # mean ocean depth
td = 180.*365.              # diffusion time scale
ta = 29.*365.               # advective time scale
q = 192e10*3600.*24.        # transport coefficient
V = 300*4.5*8200.*1e9       # Ocean volume
alphat = 1e-4               # thermal expansion coefficient
alphas = 76e-5              # haline contraction coefficient
S0 = 35                     # Reference salinity
Th = 25                     # meridional temp diff
rho0 = 1029.                # Reference density
F0 = 2.3/365                # Scale freshwater flux
        
# ---------------------------------------------- #
# INITIAL CONDITIONS
# ---------------------------------------------- #
        
DT0 = 0.2
DS0 = 0.2

# ---------------------------------------------- #
# HELPER EQUATIONS
# ---------------------------------------------- #

def Q(Drho):
    return 1./td+q*Drho**2./(rho0**2*V)
    
def Drhof(DT,DS):
    return 1.-alphat*DT+alphas*DS
    
def Fs(t):
    return min(F0+max(t-200,0)*0.1*F0,5*F0)

# ---------------------------------------------- #
# DIFFERENTIAL EQUATIONS
# ---------------------------------------------- #

def dDT(t,DT,DS):
    return -1./tr*(DT-Th)-Q(Drhof(DT,DS))*DT

def dDS(t,DT,DS):
    return Fs(t)*S0/H-Q(Drhof(DT,DS))*DS

# ---------------------------------------------- #
# TIME INTEGRATION LOOP
# ---------------------------------------------- #

DTvec=[DT0]
DSvec=[DS0]
Fvec=[Fs(0)]
tvec=[0]
i=0
while i<tmax/dt:
    t=i*dt
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
    
    Fvec.append(Fs(t))
    tvec.append(t)
    i=i+1

tvec=np.array(tvec)/365.

plt.plot(tvec,DTvec,linewidth=3)
plt.plot(tvec,DSvec,linewidth=3)
plt.plot(tvec,np.array(Fvec)*10,linewidth=3)
plt.legend([r'$\Delta T$',r'$\Delta S$',r'$Fs(t)$'],loc='best')
plt.xlim([0,tvec[-1]])
plt.ylim([0,30])