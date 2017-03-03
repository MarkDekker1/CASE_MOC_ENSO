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

tmax=1000000.
dt=10

# ---------------------------------------------- #
# PHYSICAL PARAMETERS (TAB 10.3, NCD) (IN DAYS)
# ---------------------------------------------- #

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
        
# ---------------------------------------------- #
# INITIAL CONDITIONS
# ---------------------------------------------- #
        
DT0 = 0.2
DS0 = 0.2

# ---------------------------------------------- #
# FORCING
# ---------------------------------------------- #

def force():
    global Fs
    Fs=min(0.006+max(0.00000005*(t-400000),0),0.0082)
    
# ---------------------------------------------- #
# HELPER EQUATIONS
# ---------------------------------------------- #
    
def flow(DT,DS):
    return gamma*Drhof(DT,DS)/rho0

def Q(Drho):
    return 1./td+q*Drho**2./(V)
    
def Drhof(DT,DS):
    return rho0*(alphat*DT-alphas*DS)
    
# ---------------------------------------------- #
# DIFFERENTIAL EQUATIONS
# ---------------------------------------------- #

def dDT(t,DT,DS):
    return -1./tr*(DT-Th)-Q(Drhof(DT,DS))*DT

def dDS(t,DT,DS):
    return Fs*S0/H-Q(Drhof(DT,DS))*DS

# ---------------------------------------------- #
# TIME INTEGRATION LOOP
# ---------------------------------------------- #

DTvec=[DT0]
DSvec=[DS0]
PSIvec=[DS0]
Fvec=[Fs]
PSIvec=[flow(DT0,DS0)]
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
    
    Fvec.append(Fs)
    PSIvec.append(flow(DTnew,DSnew))
    tvec.append(t)
    force()
    i=i+1

tvec=np.array(tvec)/1.
#plt.plot(tvec,DTvec,linewidth=3)
#plt.plot(tvec,DSvec,linewidth=3)
plt.plot(tvec,np.array(Fvec),linewidth=3)
plt.plot(tvec,np.array(PSIvec),linewidth=3)
plt.legend([r'$F_s$',r'$\Psi$',r'$Fs(t)*100$'],loc='best')
plt.xlim([10000,tvec[-1]])
#plt.ylim([-0.002,0.02])
plt.xlabel('Time',fontsize=15)
plt.ylabel(r'Variables',fontsize=15)
plt.tick_params(axis='both',which='major',labelsize=15)