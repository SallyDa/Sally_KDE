
import numpy as np
import astropy.units as u
import copy

def rad(submap):
    '''
    a function that calculates the heliocentric earth equatorial coordinates phi & theta
    '''
    r_sun = submap.meta['rsun_obs']*u.arcsec
    x = np.arange(submap.dimensions[0].value)
    y = np.arange(submap.dimensions[1].value)
    xx, yy = np.meshgrid(x, y)
    sx, sy = submap.pixel_to_data(xx*u.pix, yy*u.pix)
    phi = np.arcsin(sx/r_sun)
    theta = np.arcsin(sy/(r_sun*np.cos(phi)))

    phi = np.nan_to_num(phi/u.rad) * u.rad # replace nan (outside the Sun) in the array by zero
    theta = np.nan_to_num(theta/u.rad) * u.rad

    return [phi.to(u.deg), theta.to(u.deg)]


def map_rad(submap):
    '''
    a function that applies the cosine correction to the submap data
    returns a radialised map
    '''
    phi, theta = rad(submap)
    rmap = copy.deepcopy(submap)
    rmap.data = rmap.data/(np.cos(phi)*np.cos(theta))

    return rmap
    