
import os.path
import tempfile

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sunpy.map
import astropy.units as u

from Sally_KDE import radialise, kde, slope

data_directory = tempfile.mkdtemp(prefix='KDE_example')

### download data ###
from sunpy.net import vso
client = vso.VSOClient()
result = client.query(vso.attrs.Time((2015, 11, 3, 13, 30), (2015, 11, 3, 13, 31)), 
                      vso.attrs.Instrument('hmi'), vso.attrs.Physobs('LOS_magnetic_field'))
download = client.get(result, path=os.path.join(data_directory, 'decaying_region.fits')).wait()
# Note: to save the file with its proper file name, use path=os.path.join(data_directory, '{file}')

### open data file, make submap and radialse ###
mymap = sunpy.map.Map(os.path.join(data_directory, 'decaying_region.fits')).rotate()
smap = mymap.submap(u.quantity.Quantity([-250, 250]*u.arcsec),
                    u.quantity.Quantity([-450, -200]*u.arcsec))
rsmap = radialise.map_rad(smap) # radialise submap

### plot figure ###
fig = plt.figure(figsize=(6, 10))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.25])
ax1 = plt.subplot(gs[0])
rsmap.plot(vmin=-200, vmax=200) # plot the submap
ax2 = plt.subplot(gs[1])
kde1 = kde.kde_plot(rsmap, ax=ax2) # plot the kde
# to find the knees, calculate first the extra smooth kde
kde2data = kde.kde_plot(rsmap, ax=None, extrasmooth=2)[0]
slope2 = slope.slope_plot(kde2data) # find the slope of the kde plot
# find the knee positions, based on the extra smoothed kde data
knee1p = slope.B_at_steepest_slope(slope.cut_slope_array(slope2[1], cutoff=-1))
knee1n = slope.B_at_steepest_slope(slope.cut_slope_array(slope2[2], cutoff=-1))
knee2p = slope.B_at_second_knee(slope2[1], cutoff=-3)
knee2n = slope.B_at_second_knee(slope2[2], cutoff=-3)
# use the data from the original kde (not extra smooth) to find the best fit lines 
# split the data array from large negative to large positive into two absolute arrays
parray, narray = slope.split_kde_data_array(kde1[1])
# go in from the knee positions and find best fit lines between these points
Bminp = knee1p * 1.5
Bminn = knee1n * 1.5
Bmaxp = knee2p / 1.5
Bmaxn = knee2n / 1.5
parray_bknees = slope.cut_array(parray, Bminp, Bmaxp)
bfl_p = slope.bfl_loglog(parray_bknees, ax=ax2, color='r')
slopep = bfl_p[1]
narray_bknees = slope.cut_array(narray, Bminn, Bmaxn)
bfl_n = slope.bfl_loglog(narray_bknees, ax=ax2, color='b')
slopen = bfl_n[1]
print('positive slope {}, negative slope {}'.format(slopep, slopen))
# plot vertical lines to indicate the position of the knees
fvals1 = np.asarray([10**-6, 10**-1])
for kneep in [knee1p, knee2p]:
    v_pos = np.array([1, 1])*kneep
    ax2.plot(v_pos, fvals1, color='r', alpha=0.3)
for kneen in [knee1n, knee2n]:
    v_neg = np.array([1, 1])*kneen
    ax2.plot(v_neg, fvals1, color='b', alpha=0.3)
plt.tight_layout()  
plt.show()
