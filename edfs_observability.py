import matplotlib.pyplot as plt

import numpy as np
import scipy.ndimage
import scipy.interpolate
import pandas as pd

from astropy import units as u
import astropy.time
import astropy.coordinates

plt.ion()

##########

class Target:
    
    def __init__(self, ra, dec, name='None', area=None):
        """
        Specify target in equatorial coordinates in decimal degrees
        """
        self.coord = astropy.coordinates.SkyCoord(float(ra)*u.deg, 
                                                  float(dec)*u.deg, 
                                                  frame='icrs')
        self.name = name
        self.area = area
        
    def getAirmass(self, time, location):
        """
        Compute airmass of target for a given time and location
        """
        altaz_frame = astropy.coordinates.AltAz(obstime=time, location=location)
        airmass = self.coord.transform_to(altaz_frame).secz.value
        airmass[airmass < 1.] = np.inf
        return airmass
    
    def computeAirmassArray(self, time_array, location):
        """
        Compute airmas
        """
        self.time_array = time_array
        self.airmass_array = self.getAirmass(time_array, location)

    #def getOptimalTime(self, time_array, location):
    #    self.computeAirmassArray(time_array, location)
    #    return self.time_array[np.argmin(self.airmass_array)]

##########
        
def plotAirmass(time_array, airmass_array, **kwargs):
    time_array = np.array(time_array)
    airmass_array = np.array(airmass_array)
    
    label_array, n_labels = scipy.ndimage.label(airmass_array > 1.)
    for label in range(0, n_labels):
        cut = (label_array == label + 1)
        plt.plot(time_array[cut], airmass_array[cut], **kwargs)
        
##########
        
def getSunAltitude(time, location):
    altaz_frame = astropy.coordinates.AltAz(obstime=time, location=location)
    altitude = astropy.coordinates.get_sun(time).transform_to(altaz_frame).alt.degree
    return altitude

##########

def getMoonAltitude(time, location):
    altaz_frame = astropy.coordinates.AltAz(obstime=time, location=location)
    altitude = astropy.coordinates.get_moon(time).transform_to(altaz_frame).alt.degree
    return altitude

##########

# https://en.wikipedia.org/wiki/Vera_C._Rubin_Observatory
#cerro_pachon = astropy.coordinates.EarthLocation(lat=-30.244639*u.deg, lon=-70.749417*u.deg, height=2663*u.m)

# http://www.ctio.noao.edu/noao/content/coordinates-observatories-cerro-tololo-and-cerro-pachon
LON_CTIO = '-70d48m23.49s'
LAT_CTIO = '-30d10m10.78s'
ELEVATION_CTIO = 2206.8 # m
cerro_tololo = astropy.coordinates.EarthLocation(lat=LAT_CTIO,
                                                 lon=LON_CTIO,
                                                 height=ELEVATION_CTIO*u.m)

TARGET_EDFS = False
if TARGET_EDFS:
    target_array = [Target(61.24, -48.42)]
else:
    target_array = [Target(54.2743, -27.1116, 'C1'),
                    Target(54.2743, -29.0884, 'C2'),
                    Target(52.6484, -28.1000, 'C3'),
                    Target(34.4757, -4.9295, 'X1'),
                    Target(35.6645, -6.4121, 'X2'),
                    Target(36.4500, -4.6000, 'X3'),
                    Target(42.8200, 0.0000, 'S1'),
                    Target(41.1944, -0.9884, 'S2'),
                    Target(7.8744, -43.0096, 'E1'),
                    Target(9.5000, -43.9980, 'E2')]

time_search = astropy.time.Time('2020-10-17 00:00:00', scale='utc')
#time_search = astropy.time.Time('2020-10-15 00:00:00', scale='utc')
#time_search = astropy.time.Time('2020-10-14 00:00:00', scale='utc')
#time_search = astropy.time.Time('2020-10-13 00:00:00', scale='utc')
#time_search = astropy.time.Time('2020-10-29 00:00:00', scale='utc')
time_array = time_search + np.linspace(-12., 12., (16 * 24) + 1) * u.hr

for ii, target in enumerate(target_array):
    target.computeAirmassArray(time_array, cerro_tololo)

tick_values = np.linspace(np.min(time_array.decimalyear),
                          np.max(time_array.decimalyear),
                          24 + 1)
tick_labels = []
for tick_value in tick_values:
    time = astropy.time.Time(tick_value, format='decimalyear', scale='utc')
    #tick_labels.append(time.iso.split()[1].split('.')[0])
    tick_labels.append(':'.join(time.iso.split()[1].split(':')[0:2]))

sun_altitude = getSunAltitude(time_array, cerro_tololo)
moon_altitude = getMoonAltitude(time_array, cerro_tololo)

plt.figure(dpi=100)
plt.fill_between(target.time_array.decimalyear, 1., 2., where=(sun_altitude > -14.), zorder=0, color='red', alpha=0.1)
plt.fill_between(target.time_array.decimalyear, 1., 2., where=((sun_altitude < -14.) & (moon_altitude < 0.)), zorder=0, color='green', alpha=0.1)
plt.fill_between(target.time_array.decimalyear, 1., 2., where=((sun_altitude < -14.) & (moon_altitude > 0.)), zorder=0, color='yellow', alpha=0.1)
plt.axvline(time_search.decimalyear, c='black', ls='-', lw=1)
for target in target_array:
    plotAirmass(target.time_array.decimalyear, target.airmass_array, c='black', lw=1, alpha=0.25)
    index = np.argmin(target.airmass_array)
    optimal_airmass = target.airmass_array[index]
    optimal_time = target.time_array[index].decimalyear
    if '1' in target.name:
        plt.text(optimal_time, optimal_airmass + 0.02, target.name,
                 horizontalalignment='center', verticalalignment='center')
    print('%10s: %s'%(target.name, target.time_array[index]))
plt.xlim(np.min(time_array.decimalyear),
         np.max(time_array.decimalyear))
plt.ylim(1., 2.)
plt.ylabel('Airmass')
plt.xlabel('Time (UTC)')
plt.xticks(tick_values, tick_labels)
plt.title(time_search.iso)

print(target.time_array[sun_altitude < -14.][0])
print(target.time_array[sun_altitude < -14.][-1])

##########

# Check lunar distance and illumination

import ephem

#lunar_separation = astropy.coordinates.get_moon(time_array).separation(target.coord)
#print(np.min(lunar_separation))

#target_array = [Target(61.24, -48.42)]

#time_search = astropy.time.Time('2020-10-13 00:00:00', scale='utc')
time_array = time_search + np.arange(-12., 12. + (24. * 6), 0.1) * u.hr

tick_values = np.linspace(np.min(time_array.decimalyear),
                          np.max(time_array.decimalyear),
                          5)
tick_labels = []
for tick_value in tick_values:
    time = astropy.time.Time(tick_value, format='decimalyear', scale='utc')
    tick_labels.append(time.iso.split()[1].split('.')[0])

sun_altitude = getSunAltitude(time_array, cerro_tololo)
moon_altitude = getMoonAltitude(time_array, cerro_tololo)
lunar_separation = astropy.coordinates.get_moon(time_array).separation(target.coord)

moon_phase = np.empty(len(time_array))
for ii in range(0, len(time_array)):
    moon = ephem.Moon(time_array[ii].datetime)
    moon_phase[ii] = moon.moon_phase

plt.figure(dpi=100)
plt.fill_between(time_array.decimalyear, 0., 180., where=(sun_altitude > -14.), zorder=0, color='red', alpha=0.1)
plt.fill_between(time_array.decimalyear, 0., 180., where=((sun_altitude < -14.) & (moon_altitude < 0.)), zorder=0, color='green', alpha=0.1)
plt.fill_between(time_array.decimalyear, 0., 180., where=((sun_altitude < -14.) & (moon_altitude > 0.)), zorder=0, color='yellow', alpha=0.1)
plt.axvline(time_search.decimalyear, c='black', ls='-', lw=1)
#for target in target_array:
#    plotAirmass(target.time_array.decimalyear, target.airmass_array, c='black', lw=1, alpha=0.25)
plt.scatter(time_array.decimalyear, lunar_separation, c=moon_phase, vmin=0, vmax=1, cmap='binary_r')
plt.colorbar(label='Lunar Illumination')
#plt.plot(time_array.decimalyear, moon_phase, c='black')
plt.xlim(np.min(time_array.decimalyear),
         np.max(time_array.decimalyear))
plt.ylim(0., 180.)
plt.ylabel('Lunar Separation (deg)')
plt.xlabel('Time (UTC)')
plt.xticks(tick_values, tick_labels)
plt.title(time_search.iso)



"""
plt.figure(dpi=100)
plt.fill_between(target.time_array.decimalyear, 1., 2., where=(sun_altitude > -14.), zorder=0, color='red', alpha=0.1)
plt.fill_between(target.time_array.decimalyear, 1., 2., where=((sun_altitude < -14.) & (moon_altitude < 0.)), zorder=0, color='green', alpha=0.1)
plt.fill_between(target.time_array.decimalyear, 1., 2., where=((sun_altitude < -14.) & (moon_altitude > 0.)), zorder=0, color='yellow', alpha=0.1)
plt.axvline(time_search.decimalyear, c='black', ls='-', lw=1)
for target in target_array:
    plotAirmass(target.time_array.decimalyear, target.airmass_array, c='black', lw=1, alpha=0.25)
plt.xlim(np.min(time_array.decimalyear),
         np.max(time_array.decimalyear))
plt.ylim(1., 2.)
plt.ylabel('Airmass')
plt.xlabel('Time (UTC)')
plt.xticks(tick_values, tick_labels)
plt.title(time_search.iso)
"""

##########
"""
import ephem
moon = ephem.Moon(time_search.datetime)
moon.moon_phase
#optimal_time = 
#moon = ephem.Moon(optimal_time)
"""
"""
print 'Event'
print '  Event ID = %s'%(self.eventid)
print '  (ra, dec) = (%.2f, %.2f)'%(self.ra, self.dec)

print 'Date'
print '  Now = %s (UTC)'%(ephem.now().__str__())
print '  Search time = %s (UTC)'%(self.later.__str__())
print '  Optimal time = %s (UTC)'%(self.optimal_time.__str__())
print '  Airmass at optimal time = %.2f'%(self.minimum_airmass)

print 'Sun'
sun = ephem.Sun(self.optimal_time)
ra_sun, dec_sun = np.degrees(sun.ra), np.degrees(sun.dec)
print '  Angular separation = %.2f (deg)'%(utils.angsep(ra_target, dec_target, ra_sun, dec_sun))    
print '  Next rising = %s (UTC)'%(self.observatory.next_rising(ephem.Sun()).__str__())
print '  Next setting = %s (UTC)'%(self.observatory.next_setting(ephem.Sun()).__str__())

print 'Moon'
moon = ephem.Moon(self.optimal_time)
ra_moon, dec_moon = np.degrees(moon.ra), np.degrees(moon.dec)
print '  Illumination = %.2f'%(moon.moon_phase)
print '  Angular separation = %.2f (deg)'%(utils.angsep(ra_target, dec_target, ra_moon, dec_moon))
print '  Next rising = %s (UTC)'%(self.observatory.next_rising(ephem.Moon()).__str__())
print '  Next setting = %s (UTC)'%(self.observatory.next_setting(ephem.Moon()).__str__())
print '  Next new moon = %s (UTC)'%(ephem.next_new_moon(ephem.now()).__str__())
print '  Next full moon = %s (UTC)'%(ephem.next_full_moon(ephem.now()).__str__())
"""
