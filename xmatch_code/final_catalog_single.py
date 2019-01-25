'''
Reads in each of the catalogs from LOBES fields around a single EoR field. Cross-matches catalogs from each field to determine common sources. Reads out final catalog that contains all sources from all fields without repeating any sources.
'''

import numpy as np
import astropy.wcs as wcs
import astropy.io.fits as pf
from astropy import units as u
from astropy.coordinates import ICRS,SkyCoord
from astropy.coordinates import Angle
from astropy.table import Table
import glob
import os
import sys
from math import *

#### Functions to determine angluar distance between two points:


def great_circle_dist1(r1, d1, r2, d2):
    return acos(sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(r1 - r2))

def deg2rad(x):
    return x * np.pi / 180.0

def precise_dist(ra1, dec1, ra2, dec2):
    d1 = deg2rad(dec1)
    d2 = deg2rad(dec2)
    r1 = deg2rad(ra1)
    r2 = deg2rad(ra2)
    dist = great_circle_dist1(r1, d1, r2, d2)
    dist = dist * 180.0/np.pi
    return dist

#########################################################################################
## MAIN PROGRAM
#########################################################################################
path = os.getcwd()
chn = sys.argv[1]
fov = np.round((234.2/float(chn))*(1./4.509)*(180./np.pi), decimals=0) #deg
tbl_name = 'LOBES_chn%s.fits' %(chn)

## Read in aegean catalogs
tables = glob.glob('LOBES*_%s_comp.fits' %(chn))
tables.sort()
tables = np.array(tables)

## Create an array that contains all the relevant information from each 
## aegean catalog. Also only select sources within the image FOV
src_names = []
radegs = []
err_ra = []
decdegs = []
err_dec = []
flx= []
err_flx = []
major = []
err_major = []
minor = []
err_minor = []
pa = []
err_pa = []

for t in tables:

	data,header= pf.getdata('LOBES%s_%s.fits' %(t[5], chn), header=True)
	racen = float(header['CRVAL1']) #deg
	deccen = float(header['CRVAL2']) #deg
	
	print racen/15.
	print deccen
'''
	num_srcs = len(tbl['ra'][:])
	
	
	radegs.append(np.array(tbl['ra'][:]))
	ras = np.array(tbl['ra'][:])
	err_ra.append(np.array(tbl['err_ra'][:]))
	decdegs.append(np.array(tbl['dec'][:]))
	decs = np.array(tbl['dec'][:])
	err_dec.append(np.array(tbl['err_dec'][:]))
	flx.append(np.array(tbl['int_flux'][:]))
	err_flx.append(np.array(tbl['err_int_flux'][:]))
	major.append(np.array(tbl['a'][:]))
	err_major.append(np.array(tbl['err_a'][:]))
	minor.append(np.array(tbl['b'][:]))
	err_minor.append(np.array(tbl['err_b'][:]))
	pa.append(np.array(tbl['pa'][:]))
	err_pa.append(np.array(tbl['err_pa'][:]))
	sc_nme = []
	for n in range(num_srcs):
		text = 'aegean_'+str(np.around(ras[n], decimals=5))+'_'+str(np.around(decs[n], decimals=5))
		sc_nme.append(text)
	src_names.append(sc_nme)


lobes6 = SkyCoord(radegs[0][:]*u.degree, decdegs[0][:]*u.degree, frame='icrs')
lobes7 = SkyCoord(radegs[1][:]*u.degree, decdegs[1][:]*u.degree, frame='icrs')



## 7 & 6
indxlobes7, dist2lobes7, dist3lobes7 = lobes6.match_to_catalog_sky(lobes7)
indx7_match = []
i = 0
for d,g in enumerate(indxlobes7):
	aa = dist2lobes7[d]*u.degree
	if aa.value < 0.0272:
		indx7_match.append(g)
		i += 1

ra7 = np.delete(radegs[1][:], indx7_match)
sname7 = np.delete(src_names[1][:], indx7_match)
era7 = np.delete(err_ra[1][:], indx7_match)
dec7 = np.delete(decdegs[1][:], indx7_match)
edec7 = np.delete(err_dec[1][:], indx7_match)
flx7 = np.delete(flx[1][:], indx7_match)
eflx7 = np.delete(err_flx[1][:], indx7_match)
mjr7 = np.delete(major[1][:], indx7_match)
emjr7 = np.delete(err_major[1][:], indx7_match)
mnr7 = np.delete(minor[1][:], indx7_match)
emnr7 = np.delete(err_minor[1][:], indx7_match)
pa7 = np.delete(pa[1][:], indx7_match)
epa7 = np.delete(err_pa[1][:], indx7_match)

## 7 & 8

lobes7 = SkyCoord(ra7*u.degree, dec7*u.degree, frame='icrs')
lobes8 = SkyCoord(radegs[2][:]*u.degree, decdegs[2][:]*u.degree, frame='icrs')
indxlobes78, dist2lobes78, dist3lobes78 = lobes8.match_to_catalog_sky(lobes7)

indx78_match = []
i = 0
for d,g in enumerate(indxlobes78):
	aa = dist2lobes78[d]*u.degree
	if aa.value < 0.0272:
		indx78_match.append(g)
		i += 1

ra78 = np.delete(ra7, indx78_match)
sname78 = np.delete(sname7, indx78_match)
era78 = np.delete(era7, indx78_match)
dec78 = np.delete(dec7, indx78_match)
edec78 = np.delete(edec7, indx78_match)
flx78 = np.delete(flx7, indx78_match)
eflx78 = np.delete(eflx7, indx78_match)
mjr78 = np.delete(mjr7, indx78_match)
emjr78 = np.delete(emjr7, indx78_match)
mnr78 = np.delete(mnr7, indx78_match)
emnr78 = np.delete(emnr7, indx78_match)
pa78 = np.delete(pa7, indx78_match)
epa78 = np.delete(epa7, indx78_match)

## 6 & 8
indxlobes8, dist2lobes8, dist3lobes8 = lobes6.match_to_catalog_sky(lobes8)
indx8_match = []
i = 0
for d,g in enumerate(indxlobes8):
	aa = dist2lobes8[d]*u.degree
	if aa.value < 0.0272:
		indx8_match.append(g)
		i += 1

ra8 = np.delete(radegs[2][:], indx8_match)
sname8 = np.delete(src_names[2][:], indx8_match)
era8 = np.delete(err_ra[2][:], indx8_match)
dec8 = np.delete(decdegs[2][:], indx8_match)
edec8 = np.delete(err_dec[2][:], indx8_match)
flx8 = np.delete(flx[2][:], indx8_match)
eflx8 = np.delete(err_flx[2][:], indx8_match)
mjr8 = np.delete(major[2][:], indx8_match)
emjr8 = np.delete(err_major[2][:], indx8_match)
mnr8 = np.delete(minor[2][:], indx8_match)
emnr8 = np.delete(err_minor[2][:], indx8_match)
pa8 = np.delete(pa[2][:], indx8_match)
epa8 = np.delete(err_pa[2][:], indx8_match)



## 7 & 9
lobes7 = SkyCoord(ra78*u.degree, dec78*u.degree, frame='icrs')
lobes9 = SkyCoord(radegs[3][:]*u.degree, decdegs[3][:]*u.degree, frame='icrs')
indxlobes79, dist2lobes79, dist3lobes79 = lobes9.match_to_catalog_sky(lobes7)
indx79_match = []
i = 0
for d,g in enumerate(indxlobes79):
	aa = dist2lobes79[d]*u.degree
	if aa.value < 0.0272:
		indx79_match.append(g)
		i += 1

ra79 = np.delete(ra78, indx79_match)
sname79 = np.delete(sname78, indx79_match)
era79 = np.delete(era78, indx79_match)
dec79 = np.delete(dec78, indx79_match)
edec79 = np.delete(edec78, indx79_match)
flx79 = np.delete(flx78, indx79_match)
eflx79 = np.delete(eflx78, indx79_match)
mjr79 = np.delete(mjr78, indx79_match)
emjr79 = np.delete(emjr78, indx79_match)
mnr79 = np.delete(mnr78, indx79_match)
emnr79 = np.delete(emnr78, indx79_match)
pa79 = np.delete(pa78, indx79_match)
epa79 = np.delete(epa78, indx79_match)


## 8 & 9
lobes8 = SkyCoord(ra8*u.degree, dec8*u.degree, frame='icrs')
lobes9 = SkyCoord(radegs[3][:]*u.degree, decdegs[3][:]*u.degree, frame='icrs')
indxlobes89, dist2lobes89, dist3lobes89 = lobes9.match_to_catalog_sky(lobes8)
indx89_match = []
i = 0
for d,g in enumerate(indxlobes89):
	aa = dist2lobes89[d]*u.degree
	if aa.value < 0.0272:
		indx89_match.append(g)
		i += 1

ra89 = np.delete(ra8, indx89_match)
sname89 = np.delete(sname8, indx89_match)
era89 = np.delete(era8, indx89_match)
dec89 = np.delete(dec8, indx89_match)
edec89 = np.delete(edec8, indx89_match)
flx89 = np.delete(flx8, indx89_match)
eflx89 = np.delete(eflx8, indx89_match)
mjr89 = np.delete(mjr8, indx89_match)
emjr89 = np.delete(emjr8, indx89_match)
mnr89 = np.delete(mnr8, indx89_match)
emnr89 = np.delete(emnr8, indx89_match)
pa89 = np.delete(pa8, indx89_match)
epa89 = np.delete(epa8, indx89_match)

## 9 & 6

indxlobes9, dist2lobes9, dist3lobes9 = lobes6.match_to_catalog_sky(lobes9)
indx9_match = []
i = 0
for d,g in enumerate(indxlobes9):
	aa = dist2lobes9[d]*u.degree
	if aa.value < 0.0272:
		indx9_match.append(g)
		i += 1

ra9 = np.delete(radegs[3][:], indx9_match)
sname9 = np.delete(src_names[3][:], indx9_match)
era9 = np.delete(err_ra[3][:], indx9_match)
dec9 = np.delete(decdegs[3][:], indx9_match)
edec9 = np.delete(err_dec[3][:], indx9_match)
flx9 = np.delete(flx[3][:], indx9_match)
eflx9 = np.delete(err_flx[3][:], indx9_match)
mjr9 = np.delete(major[3][:], indx9_match)
emjr9 = np.delete(err_major[3][:], indx9_match)
mnr9 = np.delete(minor[3][:], indx9_match)
emnr9 = np.delete(err_minor[3][:], indx9_match)
pa9 = np.delete(pa[3][:], indx9_match)
epa9 = np.delete(err_pa[3][:], indx9_match)


srcN = np.append(src_names[0][:], np.append(sname79, np.append(sname89, sname9)))
ras = np.append(radegs[0][:], np.append(ra79, np.append(ra89, ra9)))
err_ras = np.append(err_ra[0][:], np.append(era79, np.append(era89, era9)))
decs = np.append(decdegs[0][:], np.append(dec79, np.append(dec89, dec9)))
err_decs = np.append(err_dec[0][:], np.append(edec79, np.append(edec89, edec9)))
flxs = np.append(flx[0][:], np.append(flx79, np.append(flx89, flx9)))
err_flxs = np.append(err_flx[0][:], np.append(eflx79, np.append(eflx89, eflx9)))
majors = np.append(major[0][:], np.append(mjr79, np.append(mjr89, mjr9)))
err_majors = np.append(err_major[0][:], np.append(emjr79, np.append(emjr89, emjr9)))
minors = np.append(minor[0][:], np.append(mnr79, np.append(mnr89, mnr9)))
err_minors = np.append(err_minor[0][:], np.append(emnr79, np.append(emnr89, emnr9)))
pas = np.append(pa[0][:], np.append(pa79, np.append(pa89, pa9)))
err_pas = np.append(err_pa[0][:], np.append(epa79, np.append(epa89, epa9)))

print len(ras)

tbl_out = Table([srcN, ras, err_ras, decs, err_decs, flxs, err_flxs, majors, err_majors, minors, err_minors, pas, err_pas], names=('NAMES', 'RA', 'ERR_RA', 'DEC', 'ERR_DEC', 'INT_FLUX', 'ERR_INT_FLUX', 'a', 'ERR_a', 'b', 'ERR_b', 'PA', 'ERR_PA'))
tbl_out.write(tbl_name, format='fits')
'''

