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
file = '/home/mark/Documents/LOCAL/Data/Thesis/CASEH/Data.nc'
ncdf = Dataset(file, mode='r')
Time = ncdf.variables['time'][:]
Lat  = ncdf.variables['latitude'][:]
Lon  = ncdf.variables['longitude'][:]
U    = ncdf.variables['u10'][:][:]
T    = ncdf.variables['sst'][:][:]

#Lat  = Lat[::-1]
#Lon  = list(Lon[180:])+list(Lon[0:180])
Lat  = np.array(Lat)
Lon  = np.array(Lon)
T    = np.array(T)
tau  = np.array(U)**2.*1.22*0.0013#*np.sign(U)
U    = np.array(U)

# ---------------------------------------------- #
# CLEAN VARIABLES
# ---------------------------------------------- #

a=np.where(T<265)
T[a]='nan'
tau[a]='nan'

# ---------------------------------------------- #
# SPECIFY VARIABLES
# ---------------------------------------------- #

P1 = np.where(Lat==60)[0][0]
P2 = np.where(Lat==90)[0][0]
E1 = np.where(Lat==-20)[0][0]
E2 = np.where(Lat==20)[0][0]
O1 = np.where(Lon==120)[0][0]
O2 = np.where(Lon==179)[0][0]
O3 = np.where(Lon==-140)[0][0]
O4 = np.where(Lon==-80)[0][0]

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
    #T1.append(np.nanmean(list(T[t,E2:E1,-(len(Lon)-O1):-1][0])+list(T[t,E2:E1,0:O2][0])))
    #T2.append(np.nanmean(T[t,E2:E1,O2:O3]))
    #Te.append(np.nanmean(list(T[t,E2:E1,-(len(Lon)-O1):][0])+list(T[t,E2:E1,0:O3][0])))
    #Tp.append(np.nanmean(list(T[t,P2:P1,-(len(Lon)-O1):][0])+list(T[t,P2:P1,0:O3][0])))
    T1.append(np.nanmean(T[t,E2:E1,O1:O2]))
    T2.append(np.nanmean(T[t,E2:E1,O3:O4]))
    Te.append(np.nanmean(list(T[t,E2:E1,O1:O2][0])+list(T[t,E2:E1,-180:O4][0])))
    Tp.append(np.nanmean(list(T[t,P2:P1,O1:O2][0])+list(T[t,P2:P1,-180:O4][0])))
    tauvec.append(np.nanmean(list(tau[t,E2:E1,-(len(Lon)-O1):-1][0])+list(tau[t,E2:E1,0:O3][0])))
    tauveca.append(np.nanmean(tau[t,E2:E1,:]))
    Tea.append(np.nanmean(T[t,E2:E1,:]))
    Tpa.append(np.nanmean(T[t,P2:P1,:]))
    Months.append(np.mod(t,12)+1)
    
tauvec=np.array(tauvec)
Te=np.array(Te)
Tp=np.array(Tp)
T1=np.array(T1)
T2=np.array(T2)

DT=Te-Tp
PacT=T2-T1

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



#%%
# ---------------------------------------------- #
# PLOTTING
# ---------------------------------------------- #

Lat_0=np.mean(Lat)
Lon_0=np.mean(Lon)

#Actual plot
plt.figure(num=None, figsize=(10,6),dpi=150, facecolor='w', edgecolor='k')

m = Basemap(llcrnrlon=-180,llcrnrlat=-20,urcrnrlon=180,urcrnrlat=90,
            resolution='l',projection='cyl',
            lat_ts=40,lat_0=Lat_0,lon_0=Lon_0)
												
xi, yi = m(Lon, Lat)
xi, yi = np.meshgrid(xi,yi)

Colors = m.contourf(xi,yi,tau[2],150,cmap=plt.cm.jet)

m.fillcontinents(color='grey')
m.drawparallels(np.arange(-90., 91., 30), labels=[1,0,0,0], fontsize=15)
m.drawmeridians(np.arange(-180., 180., 60), labels=[0,0,0,1], fontsize=15)
m.drawcoastlines()
cbar = m.colorbar(Colors, location='bottom', pad="30%",extend='both')
cbar.ax.tick_params(labelsize=15) 
#plt.clim([-5,5])
plt.show()

#%%
from matplotlib.mlab import PCA
data = np.array(np.random.randint(10,size=(10,3)))
data=np.transpose(np.array([DT,PacT,tauvec]))
results = PCA(data)
coefs=results.Wt
Vector1=np.dot(coefs[0],np.transpose(data))
Vector2=np.dot(coefs[1],np.transpose(data))
Vector3=np.dot(coefs[2],np.transpose(data))

print results.Wt
print 'tau correlates with PacT: '+str(round(stats.pearsonr(tauvec,PacT)[0],2))
print 'tau correlates with DT: '+str(round(stats.pearsonr(tauvec,DT)[0],2))
print 'DT correlates with PacT: '+str(round(stats.pearsonr(DT,PacT)[0],2))
#%%
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = []
y = []
z = []
for item in result.Y:
 x.append(item[0])
 y.append(item[1])
 z.append(item[2])

plt.close('all') # close all latent plotting windows
fig1 = plt.figure() # Make a plotting figure
ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
pltData = [x,y,z] 
ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
 
# make simple, bare axis lines through space:
xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
 
# label the axes 
ax.set_xlabel("x-axis label") 
ax.set_ylabel("y-axis label")
ax.set_zlabel("z-axis label")
ax.set_title("The title of the plot")
plt.show() # show the plot

#%%
from scipy.optimize import curve_fit
import scipy

def fn(x, b, c,a):
    return b*x[0] + c*x[1]+a

x = scipy.array([list(DT),list(PacT)])
y = scipy.array(tauvec)
popt, pcov = curve_fit(fn, x, y)
print popt