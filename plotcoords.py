import numpy as np
from numpy import mean, sqrt, average, square, median, std, var, sum
import astropy
from astropy.table import Table, join, vstack
import matplotlib.pyplot as pl
import matplotlib.cm as cm
from astropy.io import fits
import scipy as sp
from scipy.stats import stats
import ipdb
import pandas as pd
import glob
from glob import iglob


id=None
# number of pixels across a single chip
xpix = 2048
ypix = 4096
# number of chips in each direction
xchips = 6
ychips = 2
# number of pixels across a single bin
binsize = 128




t=Table.read("alldataf11.txt",format='ascii',names=('PTFra', 'PTFdec', 'PTFx', 'PTFy', 'PTFmag', 'PTFzeropoint', 'PTFcalmag', 'PTFmag_error',
                                 'PTF_fwhm',
                             'SDSSid', 'SDSSra', 'SDSSdec', 'SDSSmag', 'SDSStype', 'date', 'pid', 'filtnum', 'chipnum'))

xoffset=[0,2085,4168,6230,8314,10404,0,2086,4169,6253,8346,10450]
yoffset=[4135,4130,4130,4130,4130,4141,0,0,0,0,0,0,0]

deltaMag=t['PTFcalmag']-t['SDSSmag']
x0=np.ones(len(t))
y0=np.ones(len(t))
for chip in xrange(12):
    x0[t['chipnum']==chip]=xoffset[chip]
    y0[t['chipnum'] == chip] = yoffset[chip]

x=t['PTFx']+x0
y=t['PTFy']+y0

vmin=-0.1
vmax=0.1
pl.axes().set_aspect('auto', 'box')



xbins=np.arange(0,xpix*xchips,binsize)
ybins=np.arange(0,ypix*ychips,binsize)


hcount=sp.stats.binned_statistic_2d(x,y,0*deltaMag+1,statistic=np.sum,bins=[xbins,ybins])

h=sp.stats.binned_statistic_2d(x,y,deltaMag,statistic=np.nanmedian,bins=[xbins,ybins])
hs=h.statistic
#pl.imshow(hs.T,origin='lower',extent=(0,12462,0,8204), vmin=vmin, vmax=vmax)
#cbar= pl.colorbar(cmap=deltaMag, orientation='vertical')


def nanrms(b, axis=None):
    return sqrt(np.nanmean(b**2, axis=axis))

r=sp.stats.binned_statistic_2d(x,y,deltaMag,statistic=nanrms,bins=[xbins,ybins])
rs=r.statistic
pl.imshow(rs.T,origin='lower',extent=(0,12462,0,8204),vmin=0,vmax=4.6)
cbar= pl.colorbar(cmap=deltaMag, orientation='vertical')


#these are giving the medians of all bins in each "row" or "column" along either axis
#pl.plot((ybins[1:]+ybins[:-1])/2,np.nanmedian(hs,axis=0))
#pl.plot((xbins[1:]+xbins[:-1])/2,np.nanmedian(hs,axis=1))
#pl.errorbar((xbins[1:]+xbins[:-1])/2,np.nanmedian(hs,axis=1),xerr=None,yerr=std(hs,axis=1)/sqrt(len(xbins)),ecolor='black')

#these are giving the amount of x- or y-bins for a specific median
#pl.hist(np.nanmedian(hs,axis=0))
#pl.hist(np.nanmedian(hs,axis=1))


#pl.hist(nanrms(hs,axis=0))


#pl.step((ybins[1:]+ybins[:-1])/2,hs.T,where='mid')
#pl.plot((ybins[1:]+ybins[:-1])/2,np.nanmean(hs,axis=0),'o')

#these are giving me the rms of the medians of all bins in each "row" or "column"  along either axis
#pl.plot((xbins[1:]+xbins[:-1])/2, nanrms(hs,axis=1))
#pl.plot((ybins[1:]+ybins[:-1])/2, nanrms(hs,axis=0))


# this will give interquartile range
#df=pd.DataFrame(np.nanmedian(hs,axis=0))
#df.describe()

