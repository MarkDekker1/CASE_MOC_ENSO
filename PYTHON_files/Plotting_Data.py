# ---------------------------------------------- #
# PREAMBULE
# ---------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from scipy import stats
from matplotlib.mlab import PCA

#%%
# ---------------------------------------------- #
# MAP
# ---------------------------------------------- #

def boxplot(lon1,lon2,lat1,lat2,colory,style):
    x,y=m(np.zeros(lat2-lat1)+lon1,range(lat1,lat2))
    plt.plot(x,y,style,linewidth=3,color=colory)
    x,y=m(np.zeros(lat2-lat1)+lon2,range(lat1,lat2))
    plt.plot(x,y,style,linewidth=3,color=colory)
    x,y=m(range(lon1,lon2),np.zeros(lon2-lon1)+lat2)
    plt.plot(x,y,style,linewidth=3,color=colory)
    x,y=m(range(lon1,lon2),np.zeros(lon2-lon1)+lat1)
    plt.plot(x,y,style,linewidth=3,color=colory)
    

import mpl_toolkits

Lon  = ncdf.variables['longitude'][:]
velo  = ncdf.variables['si10'][:]

a=Tvar

a,Lon=mpl_toolkits.basemap.shiftgrid(0,a,Lon,start=True)

plt.figure(num=None, figsize=(10,6),dpi=150, facecolor='w', edgecolor='k')

m = Basemap(llcrnrlon=120,llcrnrlat=-10,urcrnrlon=300,urcrnrlat=90,
            resolution='l',projection='cyl',
            lat_ts=40)
												
xi, yi = m(Lon, Lat)
xi, yi = np.meshgrid(xi,yi)

Colors = m.contourf(xi,yi,a,150,cmap=plt.cm.jet,vmax=15)

m.fillcontinents(color='grey')
m.drawparallels(np.arange(-90., 91., 30), labels=[1,0,0,0], fontsize=15)
m.drawmeridians(np.arange(-180., 180., 60), labels=[0,0,0,1], fontsize=15)
m.drawcoastlines()

# Box West Pacific
boxplot(Lon1,Lon1ex,EqLat1,EqLat2,'k','-')
plt.text(Lon1*0.5+Lon1ex*0.5-3,-2,'T1',fontsize=15,color='k',weight='bold')

# Box East Pacific
boxplot(Lon2ex,Lon2,EqLat1,EqLat2,'k','-')
plt.text(Lon2ex*0.5+Lon2*0.5-3,-2,'T2',fontsize=15,color='k',weight='bold')

# Box Pole
boxplot(Lon1,Lon2,PoLat1,PoLat2,'white','--')
plt.text(0.5*Lon1+0.5*Lon2-3,80,'Tp',fontsize=15,color='white',weight='bold')

# Box Equator
boxplot(Lon1,Lon2,EqLat1,EqLat2,'white','--')
plt.text(0.5*Lon1+0.5*Lon2-3,-2,'Te',fontsize=15,color='white',weight='bold')

cbar = m.colorbar(Colors, location='bottom', pad="10%",extend='both')
cbar.ax.tick_params(labelsize=12)
plt.show()

#%%
# ---------------------------------------------- #
# Scatters
# ---------------------------------------------- #

VAR1=0.01*DTdev
VAR2=taudev-0.003*PacTdev#-0.02*PacT

xmina=np.min(VAR1)
xmaxa=np.max(VAR1)

xmin = xmina - (xmaxa-xmina)*0.1
xmax = xmaxa + (xmaxa-xmina)*0.1

ymina=np.min(VAR2)
ymaxa=np.max(VAR2)

ymin = ymina - (ymaxa-ymina)*0.1
ymax = ymaxa + (ymaxa-ymina)*0.1

fig, ax2=plt.subplots(figsize=(8,4))
ax2.scatter(VAR1, VAR2, c='r',s=50,alpha=0.7)
ax2.plot([-100,100],[0,0],'--',color='k',linewidth=1,zorder=-1)
ax2.plot([0,0],[-100,100],'--',color='k',linewidth=1,zorder=-1)

ax2.set_ylim([ymin,ymax])
ax2.set_xlim([xmin,xmax])
ax2.tick_params(axis='both',which='major',labelsize=15)
ax2.set_xlabel(r'Var1',fontsize=15)
ax2.set_ylabel(r'Var2',fontsize=15)

slope, intercept, r_value, p_value, std_err = stats.linregress(VAR1,VAR2)

ax2.text(xmina,ymaxa,'Pearson correlation is '+str(round(r_value,2)))
ax2.text(xmina,ymaxa- (ymaxa-ymina)*0.05,'Slope coefficient is '+str(round(slope,4)))

xvec=np.linspace(xmina,xmaxa,50)
fit= xvec*slope+intercept
    
ax2.plot(xvec,fit,'k',linewidth=2)

#%%
from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a+ b * x[0] + 0.002 * x[1]

xdata = [DT,PacT]
ydata = tauvec
popt1, pcov = curve_fit(func, xdata, ydata)
fitdata = func(xdata,popt1[0],popt1[1],popt1[2])

fig, (ax2,ax3)=plt.subplots(2,1,figsize=(12,4),sharex=True)

ax2.plot(Time,ydata,linewidth=2)
ax2.plot(Time,fitdata,linewidth=2)
ax2.tick_params(axis='both',which='major',labelsize=15)
ax2.set_ylabel(r'Wind stress',fontsize=15)
#ax2.set_xlabel(r'Time',fontsize=15)
ax2.legend([r'$\tau$',r'fit'])
ax2.set_xlim([Time[0],Time[-1]])
ax2.fill_between(Time, ydata, fitdata, where=ydata >= fitdata, facecolor='b', interpolate=True,alpha=0.4,color='b')
ax2.fill_between(Time, ydata, fitdata, where=ydata <= fitdata, facecolor='g', interpolate=True,alpha=0.4,color='g')

xdata = [DTdev,PacTdev]
ydata = taudev
popt2, pcov = curve_fit(func, xdata, ydata)
fitdata = func(xdata,popt2[0],popt2[1],popt2[2])

ax3.plot(Time,ydata,linewidth=2)
ax3.plot(Time,fitdata,linewidth=2)
ax3.tick_params(axis='both',which='major',labelsize=15)
ax3.set_ylabel(r'Wind stress',fontsize=15)
ax3.set_xlabel(r'Time',fontsize=15)
#ax3.legend([r'$\tau$',r'fit'])
ax3.set_xlim([Time[0],Time[-1]])
ax3.fill_between(Time, ydata, fitdata, where=ydata >= fitdata, facecolor='b', interpolate=True,alpha=0.4,color='b')
ax3.fill_between(Time, ydata, fitdata, where=ydata <= fitdata, facecolor='g', interpolate=True,alpha=0.4,color='g')
