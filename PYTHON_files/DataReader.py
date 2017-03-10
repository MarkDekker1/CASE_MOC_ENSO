# ---------------------------------------------- #
# COUPLED SYSTEM: MOC => ENSO
# ---------------------------------------------- #

# ---------------------------------------------- #
# PREAMBULE
# ---------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from scipy import stats
from matplotlib.mlab import PCA

# ---------------------------------------------- #
# READ VARIABLES
# ---------------------------------------------- #
#file = '/home/mark/Documents/LOCAL/Data/Thesis/CASEH/Data.nc'
file = 'C:\Users\Mark Dekker\Documents\Thesis\LOCAL\Data\Thesis\Data.nc'
ncdf = Dataset(file, mode='r')
Time = ncdf.variables['time'][:]
Lat  = ncdf.variables['latitude'][:]
Lon  = ncdf.variables['longitude'][:]
velo  = ncdf.variables['si10'][:]
T    = ncdf.variables['sst'][:][:]

Lat  = np.array(Lat)
Lon  = np.array(Lon)
T    = np.array(T)
tau  = -np.array(velo)**2.*1.22*0.0013

# ---------------------------------------------- #
# CLEAN VARIABLES
# ---------------------------------------------- #

a=np.where(T<265)
T[a]='nan'
tau[a]='nan'

# ---------------------------------------------- #
# SPECIFY BOXES
# ---------------------------------------------- #

EqLat1=-5
EqLat2=5
PoLat1=60
PoLat2=90

Lon1=160
Lon2=270
Lon1ex=210
Lon2ex=210

# ---------------------------------------------- #
# SPECIFY VARIABLES
# ---------------------------------------------- #

P1 = np.where(Lat==PoLat1)[0][0]
P2 = np.where(Lat==PoLat2)[0][0]
E1 = np.where(Lat==EqLat1)[0][0]
E2 = np.where(Lat==EqLat2)[0][0]
O1 = np.where(Lon==Lon1)[0][0]
O2 = np.where(Lon==179)[0][0]
O3 = np.where(Lon==-150)[0][0]
O4 = np.where(Lon==-90)[0][0]

Tp=[]
Te=[]
T1=[]
T2=[]
tauvec=[]
Months=[]
tauveca=[]
Tea=[]
Tpa=[]

for t in range(0,len(T)):
    l = list(T[t,E2:E1,O1:O2])
    p1 = [item for sublist in l for item in sublist]
    l = list(T[t,E2:E1,0:O3])
    p2 = [item for sublist in l for item in sublist]
    
    T1.append(np.nanmean(p1+p2))
        
    l = list(T[t,E2:E1,O3:O4])
    p1 = [item for sublist in l for item in sublist]
    
    T2.append(np.nanmean(p1))
        
    l = list(T[t,E2:E1,O1:O2])
    p1 = [item for sublist in l for item in sublist]
    l = list(T[t,E2:E1,0:O4])
    p2 = [item for sublist in l for item in sublist]
    
    Te.append(np.nanmean(p1+p2))
    
    l = list(T[t,P2:P1,O1:O2])
    p1 = [item for sublist in l for item in sublist]
    l = list(T[t,P2:P1,0:O4])
    p2 = [item for sublist in l for item in sublist]
    
    Tp.append(np.nanmean(p1+p2))
    
    l = list(tau[t,E2:E1,O1:O2])
    p1 = [item for sublist in l for item in sublist]
    l = list(tau[t,E2:E1,0:O4])
    p2 = [item for sublist in l for item in sublist]
    
    tauvec.append(np.nanmean(p1+p2))
    
    Months.append(np.mod(t,12)+1)
    
tauvec=np.array(tauvec)
Te=np.array(Te)
Tp=np.array(Tp)
T1=np.array(T1)
T2=np.array(T2)

DT=Te-Tp
PacT=T2-T1

# ---------------------------------------------- #
# AVERAGE AND VARIANCE T FIELD
# ---------------------------------------------- #

Tav=T[0]
for i in range(0,len(T[0])):
    for j in range(0,len(T[0,0])):
        Tav[i,j]=np.nanmean(T[:,i,j])
        
Tvar=T[1]
for i in range(0,len(T[0])):
    for j in range(0,len(T[0,0])):
        Tvar[i,j]=np.nanvar(T[:,i,j])

# ---------------------------------------------- #
# REMOVE SEASONALITY
# ---------------------------------------------- #

Dtau=[]
DTmeans=[]
Temppart=[]
Monthsing=np.linspace(1,12,12)
for i in range(0,12):
    meantje=np.nanmean(tauvec[np.arange(i,36,12)])
    Dtau.append(meantje)
    meantje=np.nanmean(Te[np.arange(i,36,12)]-Tp[np.arange(i,36,12)])
    DTmeans.append(meantje)
    meantje=np.nanmean(T2[np.arange(i,36,12)]-T1[np.arange(i,36,12)])
    Temppart.append(meantje)
    
taudev=np.zeros(len(tauvec))
DTdev=np.zeros(len(tauvec))
PacTdev=np.zeros(len(tauvec))
for i in range(0,len(tauvec)):
    taudev[i]=tauvec[i]-Dtau[Months[i]-1]
    DTdev[i]=DT[i]-DTmeans[Months[i]-1]
    PacTdev[i]=PacT[i]-Temppart[Months[i]-1]

# ---------------------------------------------- #
# CHECK STATISTICS
# ---------------------------------------------- #

data=np.transpose(np.array([DT,PacT,tauvec]))
results = PCA(data)
coefs=results.Wt
Vector1=np.dot(coefs[0],np.transpose(data))
Vector2=np.dot(coefs[1],np.transpose(data))
Vector3=np.dot(coefs[2],np.transpose(data))

datav=np.transpose(np.array([DTdev,PacTdev,taudev]))
resultsv = PCA(datav)
coefsv=resultsv.Wt
Vector1v=np.dot(coefsv[0],np.transpose(datav))
Vector2v=np.dot(coefsv[1],np.transpose(datav))
Vector3v=np.dot(coefsv[2],np.transpose(datav))

print '* --------------------------------------------------- *'
print '* ----- * ----- * ----- General ----- * ----- * ----- *'
print '* --------------------------------------------------- *'
print '* ----------------- Correlations -------------------- *'
print 'tau correlates with T2-T1: '+str(round(stats.pearsonr(tauvec,PacT)[0],2))
print 'tau correlates with DT: '+str(round(stats.pearsonr(tauvec,DT)[0],2))
print 'DT correlates with T2-T1: '+str(round(stats.pearsonr(DT,PacT)[0],2))
print '* ----------------- PCA results --------------------- *'
print 'Variances: ' + str(results.fracs)
print 'Eigenvectors (DT, T2-T1, tau):'
print results.Wt
print '* --------------------------------------------------- *'
print '* ----- * ----- * ---- Deviation ---- * ----- * ----- *'
print '* --------------------------------------------------- *'
print '* ----------------- Correlations -------------------- *'
print 'tau correlates with T2-T1: '+str(round(stats.pearsonr(taudev,PacTdev)[0],2))
print 'tau correlates with DT: '+str(round(stats.pearsonr(taudev,DTdev)[0],2))
print 'DT correlates with T2-T1: '+str(round(stats.pearsonr(DTdev,PacTdev)[0],2))
print '* ----------------- PCA results --------------------- *'
print 'Variances: ' + str(resultsv.fracs)
print 'Eigenvectors (DT, T2-T1, tau):'
print resultsv.Wt