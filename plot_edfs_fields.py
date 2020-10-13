import numpy as np
import pylab
import matplotlib.path
from matplotlib.collections import PolyCollection

from astropy import units as u
from astropy.coordinates import SkyCoord

import ugali.utils.projector

pylab.ion()

############################################################

def rotateFocalPlane(ccd_array, ra_center, dec_center, ra_field, dec_field):
    proj_center = ugali.utils.projector.Projector(ra_center, dec_center)
    proj_field = ugali.utils.projector.Projector(ra_field, dec_field)
    
    ccd_array_new = []
    for ii in range(0, len(ccd_array)):
        ra, dec = proj_field.imageToSphere(np.transpose(ccd_array[ii])[0], 
                                           np.transpose(ccd_array[ii])[1]) 
        x, y = proj_center.sphereToImage(ra, dec)
        ccd_array_new.append(zip(x, y))

    return ccd_array_new

############################################################

def plotFocalPlane(ccd_array, ra_center, dec_center, ra_field, dec_field, ax, color='red'):
    ccd_array_new = rotateFocalPlane(ccd_array, ra_center, dec_center, ra_field, dec_field)
    coll = PolyCollection(ccd_array_new, alpha=0.2, color=color, edgecolors='none')
    ax.add_collection(coll)

############################################################

#def plotPoly(ra, dec, ax, alpha=1., facecolors='none', edgecolors='black'):
def plotPoly(ra, dec, ax, **kwargs):
    x, y = proj_center.sphereToImage(ra, dec)
    #coll = PolyCollection([zip(x, y)], facecolors=facecolors, edgecolors=edgecolors)
    coll = PolyCollection([zip(x, y)], **kwargs)
    ax.add_collection(coll)

############################################################

def applyDither(ra, dec, x, y):
    proj = ugali.utils.projector.Projector(ra, dec)
    ra_dither, dec_dither = proj.imageToSphere(x, y)
    return ra_dither, dec_dither

############################################################

def plotRubinObs(ra_center, dec_center, ra_field, dec_field, ax, color='blue', **kwargs):
    proj_center = ugali.utils.projector.Projector(ra_center, dec_center)
    x_field, y_field = proj_center.sphereToImage(ra_field, dec_field)
    phi = np.linspace(0., 2. * np.pi, 361)
    radius = 1.75
    x_circle = x_field + radius * np.cos(phi)
    y_circle = y_field + radius * np.sin(phi)
    ax.plot(x_circle, y_circle, color=color, ls='--', **kwargs)
    
############################################################


# Octagon stadium envelope of EDFS
edfs_octagon = np.array([[63.25, -45.67],
                         [65.35, -46.10],
                         [66.40, -47.25],
                         [65.99, -48.72],
                         [59.25, -51.19],
                         [56.95, -50.82],
                         [55.90, -49.40],
                         [56.80, -47.99]])

# Spitzer
spitzer = [['4:12:17.743', '-45:36:40.30'],
           ['3:46:56.584', '-47:53:34.33'],
           ['3:43:53.570', '-49:13:57.66'],
           ['3:44:02.863', '-49:14:35.33'],
           ['3:43:47.506', '-49:20:39.90'],
           ['3:48:14.537', '-50:48:49.27'],
           ['3:48:34.589', '-50:47:57.65'],
           ['3:48:28.026', '-50:47:54.14'],
           ['3:57:00.298', '-51:15:14.01'],
           ['4:24:32.426', '-48:38:47.51'],
           ['4:25:58.153', '-47:23:52.12'],
           ['4:19:27.458', '-45:42:29.03']]
spitzer_radec = []
for ii in range(0, len(spitzer)):
    c = SkyCoord('%s %s'%(spitzer[ii][0].replace(':', ' '),
                          spitzer[ii][1].replace(':', ' ')), unit=(u.hourangle, u.deg))
    spitzer_radec.append([c.ra.value, c.dec.value])
spitzer_radec = np.array(spitzer_radec)

# VIRCAM
vircam = [['4:03:54.912', '-47:50:58.920'],
          ['4:06:22.992', '-48:54:59.040']]
vircam_radec = []
for ii in range(0, len(vircam)):
    c = SkyCoord('%s %s'%(vircam[ii][0].replace(':', ' '),
                          vircam[ii][1].replace(':', ' ')), unit=(u.hourangle, u.deg))
    vircam_radec.append([c.ra.value, c.dec.value])
vircam_radec = np.array(vircam_radec)

# Rubin
rubin_radec = np.array([[63.66, -47.59],
                        [58.95, -49.30]])

# DECam

def make_edfs_radec(mode='baseline'):
    ra0,dec0 = 60.9788, -47.8497
    ra1,dec1 = 61.5958, -48.9164
    if mode == 'baseline':
        dra  = [-0.1,+0.1]
        ddec = [+0.2,-0.4]
    elif mode == 'compact':
        dra  = [-0.1,+0.1]
        ddec = [+0.1,-0.3]
    elif mode == 'eight':
        dra  = [-0.1 + 1.,+0.1 + 1.]
        ddec = [+0.1 + 0.35,-0.3 + 0.35]
    decam_radec = np.array([[ra0 + dra[0], dec0 + ddec[0]],
                            [ra1 + dra[1], dec1 + ddec[1]]])
    #ra_shift, dec_shift = 2.5, 1.0
    if mode == 'baseline':
        ra_shift, dec_shift = 2.5, 0.9
    elif mode in ['compact', 'eight']:
        ra_shift, dec_shift = 2.0, 0.7
    decam_radec_2 = np.array([[ra0 + dra[0] + ra_shift, dec0 + ddec[0] + dec_shift],
                              [ra1 + dra[1] + ra_shift, dec1 + ddec[1] + dec_shift],
                              [ra0 + dra[0] - ra_shift, dec0 + ddec[0] - dec_shift],
                              [ra1 + dra[1] - ra_shift, dec1 + ddec[1] - dec_shift],
                              [ra0 + dra[0] + 2. * ra_shift, dec0 + ddec[0] + 2. * dec_shift],
                              [ra1 + dra[1] + 2. * ra_shift, dec1 + ddec[1] + 2. * dec_shift],
                              [ra0 + dra[0] - 2. * ra_shift, dec0 + ddec[0] - 2. * dec_shift],
                              [ra1 + dra[1] - 2. * ra_shift, dec1 + ddec[1] - 2. * dec_shift]])

    if mode == 'eight':
        decam_radec_2 = decam_radec_2[(0, 1, 2, 3, 6, 7),:]
    radec = np.rec.fromrecords(np.vstack([decam_radec,decam_radec_2]),names=['ra','dec'])
    return radec

decam_radec = make_edfs_radec(mode='compact')

"""
decam_radec = np.array([[60.9788 - 0.1, -47.8497 + 0.1],
                        [61.5958 + 0.1, -48.9164 - 0.3]])

ra_shift, dec_shift = 2.0, 0.7
decam_radec_2 = np.array([[60.9788 - 0.1 + ra_shift, -47.8497 + 0.1 + dec_shift],
                          [61.5958 + 0.1 + ra_shift, -48.9164 - 0.3 + dec_shift],
                          [60.9788 - 0.1 - ra_shift, -47.8497 + 0.1 - dec_shift],
                          [61.5958 + 0.1 - ra_shift, -48.9164 - 0.3 - dec_shift]])

ra_shift, dec_shift = 2.6, 0.2
decam_radec_3 = np.array([[60.9788 - 0.1 + ra_shift, -47.8497 + 0.1 - dec_shift],
                          [61.5958 + 0.1 - ra_shift, -48.9164 - 0.3 + dec_shift]])
"""

# Example
#c = SkyCoord('4 03 54.912 -47 50 58.920', unit=(u.hourangle, u.deg))
#ra_center, dec_center = c.ra.value, c.dec.value
ra_center, dec_center = 61.24, -48.42

data_ccd_corners = eval(''.join(open('ccd_corners_xy_fill.dat').readlines()))
ccd_array = []
for key in data_ccd_corners.keys():
    ccd_array.append(data_ccd_corners[key])

proj_center = ugali.utils.projector.Projector(ra_center, dec_center)

# Begin plotting

fig, ax = pylab.subplots(figsize=(8, 8))

# EDFS octagon
#for ii in range(0, len(edfs_octagon)):
#    x, y = proj_center.sphereToImage(edfs_octagon[ii][0], edfs_octagon[ii][1])
#    pylab.scatter(x, y, marker='o', color='yellow', edgecolor='black', s=200, zorder=1.e6)
#x, y = proj_center.sphereToImage(edfs_octagon[:,0], edfs_octagon[:,1])
#pylab.plot(x, y)
#pylab.Polygon(np.array([x, y]).T, ls='-', ec='black', lw=1)
plotPoly(edfs_octagon[:,0], edfs_octagon[:,1], ax, facecolors='none', edgecolors='black', label='EDFS')

#for ii in range(0, len(spitzer_radec)):
#    x, y = proj_center.sphereToImage(spitzer_radec[ii][0], spitzer_radec[ii][1])
#    pylab.scatter(x, y, marker='o', color='green', edgecolor='black', s=200, zorder=1.e6)
plotPoly(spitzer_radec[:,0], spitzer_radec[:,1], ax, facecolors='none', edgecolors='purple', label='Spitzer')

#for ii in range(0, len(vircam_radec)):
#    x, y = proj_center.sphereToImage(vircam_radec[ii][0], vircam_radec[ii][1])
#    pylab.scatter(x, y, marker='o', color='red', edgecolor='black', s=200, zorder=1.e6)
x, y = proj_center.sphereToImage(vircam_radec[:,0], vircam_radec[:,1])
pylab.scatter(x, y, marker='o', color='red', edgecolor='black', s=200, zorder=1.e6, label='VISTA/VIRCAM Centers')

#plotFocalPlane(ccd_array, ra_center, dec_center, ra_center, dec_center, ax, color='red')
for ii in range(0, len(decam_radec)):
    plotFocalPlane(ccd_array, ra_center, dec_center, decam_radec[ii][0], decam_radec[ii][1], ax, color='red')

    # Labels
    label = 'hex_%i'%(ii + 1)
    x, y = proj_center.sphereToImage(decam_radec[ii][0], decam_radec[ii][1])
    pylab.text(x, y, label, horizontalalignment='center')

    print('%s, %.2f, %.2f'%(label, decam_radec[ii][0], decam_radec[ii][1]))
    
#for ii in range(0, len(decam_radec_2)):
#    plotFocalPlane(ccd_array, ra_center, dec_center, decam_radec_2[ii][0], decam_radec_2[ii][1], ax, color='green')
#for ii in range(0, len(decam_radec_3)):
#    plotFocalPlane(ccd_array, ra_center, dec_center, decam_radec_3[ii][0], decam_radec_3[ii][1], ax, color='green')

for ii in range(0, len(rubin_radec)):
    if ii == 0:
        label = 'Rubin Obs'
    else:
        label = None
    plotRubinObs(ra_center, dec_center, rubin_radec[ii][0], rubin_radec[ii, 1], ax, color='blue', label=label)
    
pylab.xlim(4., -4.)
pylab.ylim(-4., 4.)
pylab.xlabel('$\Delta$RA (deg)', labelpad=20)
pylab.ylabel('$\Delta$Dec (deg)')
pylab.title('Centered on (RA, Dec) = (%.3f, %.3f)'%(ra_center, dec_center))

pylab.legend(loc='upper right')




##### Let's check dither #####

#decam_dither_radec = [[60.8800,  -47.7500],
#                      [60.8634, -47.7323],
#                      [60.9104, -47.7824]]
decam_dither_radec = [[60.8800, -47.7500],
                      [60.8634, -47.7323],
                      [60.8966, -47.7176]] 

fig, ax = pylab.subplots(figsize=(8, 8))

for ii in range(0, len(decam_dither_radec)):
    plotFocalPlane(ccd_array, ra_center, dec_center, decam_dither_radec[ii][0], decam_dither_radec[ii][1], ax, color='red')

pylab.xlim(4., -4.)
pylab.ylim(-4., 4.)
pylab.xlabel('$\Delta$RA (deg)', labelpad=20)
pylab.ylabel('$\Delta$Dec (deg)')
