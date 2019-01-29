'''
Reads in each of the catalogs from LOBES fields around a single EoR field. Cross-matches catalogs from each field to determine common sources. Reads out final catalog that contains all sources from all fields without repeating any sources.

Run script: python final_catalog_single.py 'chn'
'''

import numpy as np
import astropy.io.fits as pf
from astropy import units as u
from astropy.coordinates import ICRS,SkyCoord
from astropy.table import Table
import glob
import sys
from math import *

#### Functions to determine angluar distance between two points
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

#### Does cross-match between two catalogs of sources
def xmatch(rr0, dd0, rr1, dd1):
	lobes0 = SkyCoord(rr0*u.degree, dd0*u.degree, frame='icrs')
	lobes1 = SkyCoord(rr1*u.degree, dd1*u.degree, frame='icrs')
	
	indxlobes1, dist2lobes1, dist3lobes1 = lobes0.match_to_catalog_sky(lobes1)
	indx1_match = []
	i = 0
	for d,g in enumerate(indxlobes1):
		aa = dist2lobes1[d]*u.degree
		if aa.value < 0.0272:
			indx1_match.append(g)
			i += 1
	return indx1_match

## Deletes array entries associated with dindx and returns new array without entries
def del_indx(sn1, ra1, er1, de1, ede1, f1, ef1, mj1, emj1, mn1, emn1, pa1, epa1, dindx):
	
	sn_new  = np.delete(sn1  , dindx)
	ra_new  = np.delete(ra1   , dindx)
	er_new  = np.delete(er1 , dindx)
	de_new  = np.delete(de1  , dindx)
	ede_new = np.delete(ede1 , dindx)
	f_new   = np.delete(f1   , dindx)
	ef_new  = np.delete(ef1  , dindx)
	mj_new  = np.delete(mj1  , dindx)
	emj_new = np.delete(emj1 , dindx)
	mn_new  = np.delete(mn1  , dindx)
	emn_new = np.delete(emn1 , dindx)
	pa_new  = np.delete(pa1  , dindx)
	epa_new = np.delete(epa1 , dindx)
	
	return sn_new, ra_new, er_new, de_new, ede_new, f_new, ef_new, mj_new, emj_new, mn_new, emn_new, pa_new, epa_new


## Function to go through and cross match each field and return arrays with common sources removed
def xmatch_array(sn1, ra1, er1, de1, ede1, f1, ef1, mj1, emj1, mn1, emn1, pa1, epa1):
	srcnew = []
	ranew = []
	eranew = []
	decnew = []
	edecnew = []
	flxnew = []
	eflxnew = []
	mjrnew = []
	emjrnew = []
	mnrnew = []
	emnrnew = []
	panew = []
	epanew = []
	
	for i in range(len(ra1)-1):
		k = i + 1
		print 'Cross match: '+str(0)+' with '+str(k)
		indx = xmatch(ra1[0][:], de1[0][:], ra1[k][:], de1[k][:])
		snt, rat, ert, dt, edt, ft, eft, mjt, emjt, mnt, emnt, pat, ept = \
		del_indx(sn1[k][:], ra1[k][:], er1[k][:], de1[k][:], ede1[k][:], f1[k][:], ef1[k][:], mj1[k][:], emj1[k][:], mn1[k][:], emn1[k][:], pa1[k][:], epa1[k][:], indx)
		srcnew.append(snt)
		ranew.append(rat)
		eranew.append(ert)
		decnew.append(dt)
		edecnew.append(edt)
		flxnew.append(ft)
		eflxnew.append(eft)
		mjrnew.append(mjt)
		emjrnew.append(emjt)
		mnrnew.append(mnt)
		emnrnew.append(emnt)
		panew.append(pat)
		epanew.append(ept)

	return srcnew, ranew, eranew, decnew, edecnew, flxnew, eflxnew, mjrnew, emjrnew, mnrnew, emnrnew, panew, epanew



#########################################################################################
## MAIN PROGRAM
#########################################################################################
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
rad = []
errRA = []
decd = []
errDEC = []
snu= []
errSnu = []
mjr = []
errMjr = []
mnr = []
errMnr = []
posAng = []
errPosAng = []

for t in range(len(tables)):
	print tables[t]
	data,header= pf.getdata('LOBES%s_%s.fits' %(tables[t][5], chn), header=True)
	racen = float(header['CRVAL1']) #deg
	deccen = float(header['CRVAL2']) #deg
	
	radegs    = []
	err_ra    = []
	decdegs   = []
	err_dec   = []
	flx       = []
	err_flx   = []
	major     = []
	err_major = []
	minor     = []
	err_minor = []
	pa        = []
	err_pa    = []
	sc_nme    = []

	tbl = Table.read(tables[t])
	num_srcs = len(tbl['ra'][:])
	for n in range(num_srcs):
		aa1 = tbl['ra'][n]
		bb1 = tbl['dec'][n]
		dist = precise_dist(racen, deccen, aa1, bb1)
		if dist < fov:
			radegs.append(tbl['ra'][n])
			err_ra.append(tbl['err_ra'][n])
			decdegs.append(tbl['dec'][n])
			err_dec.append(tbl['err_dec'][n])
			flx.append(tbl['int_flux'][n])
			err_flx.append(tbl['err_int_flux'][n])
			major.append(tbl['a'][n])
			err_major.append(tbl['err_a'][n])
			minor.append(tbl['b'][n])
			err_minor.append(tbl['err_b'][n])
			pa.append(tbl['pa'][n])
			err_pa.append(tbl['err_pa'][n])
			text = 'aegean_'+str(np.around(tbl['ra'][n], decimals=5))+'_'+str(np.around(tbl['dec'][n], decimals=5))
			sc_nme.append(text)

	src_names.append(sc_nme)
	rad.append(radegs)
	errRA.append(err_ra)
	decd.append(decdegs)
	errDEC.append(err_dec)
	snu.append(flx)
	errSnu.append(err_flx)
	mjr.append(major)
	errMjr.append(err_major)
	mnr.append(minor)
	errMnr.append(err_minor)
	posAng.append(pa)
	errPosAng.append(err_pa)


## Start collecting final source list for each of the arrays; we will keep 
## all sources in the field 6 catalog

srcN       = src_names[0][:]
ras        = rad[0][:]
err_ras    = errRA[0][:]
decs       = decd[0][:]
err_decs   = errDEC[0][:]
flxs       = snu[0][:]
err_flxs   = errSnu[0][:]
majors     = mjr[0][:]
err_majors = errMjr[0][:]
minors     = mnr[0][:]
err_minors = errMnr[0][:]
pas        = posAng[0][:]
err_pas    = errPosAng[0][:]

## Cross-match field 6 with other fields
sc0, ra0, era0, dec0, edec0, flx0, eflx0, mjr0, emjr0, mnr0, emnr0, pa0, epa0 =\
xmatch_array(src_names, rad, errRA, decd, errDEC, snu, errSnu, mjr, errMjr, mnr, errMnr, posAng, errPosAng)

## Add all sources that are not common between field 7 and 6
srcN       = np.append(srcN, sc0[0])
ras        = np.append(ras, ra0[0])
err_ras    =  np.append(err_ras, era0[0])
decs       = np.append(decs, dec0[0])
err_decs   = np.append(err_decs, edec0[0])
flxs       = np.append(flxs, flx0[0])
err_flxs   = np.append(err_flxs, eflx0[0])
majors     = np.append(majors, mjr0[0])
err_majors = np.append(err_majors, emjr0[0])
minors     = np.append(minors, mnr0[0])
err_minors = np.append(err_minors, emnr0[0])
pas 	   = np.append(pas, pa0[0])
err_pas    = np.append(err_pas, epa0[0])


## Cross-match field 7 with other fields (8 & 9)
sc1, ra1, era1, dec1, edec1, flx1, eflx1, mjr1, emjr1, mnr1, emnr1, pa1, epa1 = \
xmatch_array(sc0, ra0, era0, dec0, edec0, flx0, eflx0, mjr0, emjr0, mnr0, emnr0, pa0, epa0)


## Add all sources that are not common between field 8 and 7 
srcN       = np.append(srcN, sc1[0])
ras        = np.append(ras, ra1[0])
err_ras    =  np.append(err_ras, era1[0])
decs       = np.append(decs, dec1[0])
err_decs   = np.append(err_decs, edec1[0])
flxs       = np.append(flxs, flx1[0])
err_flxs   = np.append(err_flxs, eflx1[0])
majors     = np.append(majors, mjr1[0])
err_majors = np.append(err_majors, emjr1[0])
minors     = np.append(minors, mnr1[0])
err_minors = np.append(err_minors, emnr1[0])
pas 	   = np.append(pas, pa1[0])
err_pas    = np.append(err_pas, epa1[0])

print len(ras)

'''
## Cross-match field 8 with field 9
sc2, ra2, era2, dec2, edec2, flx2, eflx2, mjr2, emjr2, mnr2, emnr2, pa2, epa2 = \
xmatch_array(sc1, ra1, era1, dec1, edec1, flx1, eflx1, mjr1, emjr1, mnr1, emnr1, pa1, epa1)

## Add all sources not common between field 9 and 8
srcN       = np.append(srcN, sc2)
ras        = np.append(ras, ra2)
err_ras    =  np.append(err_ras, era2)
decs       = np.append(decs, dec2)
err_decs   = np.append(err_decs, edec2)
flxs       = np.append(flxs, flx2)
err_flxs   = np.append(err_flxs, eflx2)
majors     = np.append(majors, mjr2)
err_majors = np.append(err_majors, emjr2)
minors     = np.append(minors, mnr2)
err_minors = np.append(err_minors, emnr2)
pas 	   = np.append(pas, pa2)
err_pas    = np.append(err_pas, epa2)

print 'Total number of unique sources: '+str(len(ras))

## Write out unique sources to fits table
tbl_out = Table([srcN, ras, err_ras, decs, err_decs, flxs, err_flxs, majors, err_majors, minors, err_minors, pas, err_pas], names=('NAMES', 'RA', 'ERR_RA', 'DEC', 'ERR_DEC', 'INT_FLUX', 'ERR_INT_FLUX', 'a', 'ERR_a', 'b', 'ERR_b', 'PA', 'ERR_PA'))
tbl_out.write(tbl_name, format='fits')
'''

