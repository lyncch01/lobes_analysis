import numpy as np
#import astropy.wcs as wcs
import astropy.io.fits as pf
#from astropy import units as au
#from astropy.coordinates import SkyCoord
#from astropy.coordinates import Angle
#from astropy.table import Table
import os
import glob
#import datetime
#import pytz
import sys
#from math import *
import gc
#import aplpy



################## MedianImage.py #####################################################   

gc.enable()

def cubeCreate(fileN):
	print fileN
	img = None
	gc.collect()
	img = []
	gc.collect()
	aa = 0
	for fname in fileN:
		print fname
		hdulist = pf.open(fname)
		data = np.ma.masked_array(hdulist[0].data, fill_value=np.nan)
		data[np.isnan(data)] = np.ma.masked
		#data = data[:,:,99:6499, 99:6499]
		#print np.shape(data)
		img.append(data.copy())
		hdulist.close()
		del hdulist[0].data
		data = None
		hdulist = None
		gc.collect()
        gc.collect()
	#print len(img[0])
	return img
	
def headerChanging(iname, fname):
	headerInput = pf.getheader(iname, 0)
	data, header = pf.getdata(fname, header=True)
	pf.writeto(fname, data, headerInput, clobber=True)

################## catalog create #####################################################   

def bane_run(files, c):
	if not os.path.isfile(files[:-5]+'_bkg.fits'):
		bane_cmd = 'BANE --cores=%s --compress %s' %(c, files)
		print bane_cmd
		os.system(bane_cmd)
	else:
		pass

def aegean_run(files, mj, mn, p, c, tscope, outfile):
	if not os.path.isfile(outfile+'_comp.fits'):
		aegean_cmd = 'aegean --telescope=%s --beam %s %s %s --seedclip 5 --floodclip 4 --autoload --cores=%s --table %s %s' %(tscope, mj, mn, p, c, outfile+'_comp.fits', files)
		print aegean_cmd
		os.system(aegean_cmd)
		cmd2 = 'mv '+outfile+'_comp_comp.fits '+outfile+'_comp.fits'
		os.system(cmd2)
		print cmd2
	else:
		pass


#######################################################################################
#				MAIN PROGRAM       
#######################################################################################

#set parameters of LOBES field:
chn = '093'
fld = '7'

outfile_I = 'LOBES%s_%s.fits' %(fld, chn)
fname = '%s' %(outfile_I[:-5])


#read in all files for channel
files = glob.glob('*MFS*.regrid.fits')
files.sort()

data, header = pf.getdata(files[0], header=True)
bmaj = header['BMAJ']
bmin = header['BMIN']
bpa = header['BPA']

#create a median image using all files of same channel
concat_image_I = cubeCreate(files)
print 'Created Stokes I median image: %s' %(outfile_I)

median_image_I = np.median(concat_image_I, axis=0)
hdu_I = pf.PrimaryHDU(median_image_I)
hdu_I.writeto(outfile_I, clobber=True)
headerChanging(files[0], outfile_I)

#create source catalog using aegean 
bane_run(outfile_I, '1')
print 'Run BANE'
aegean_run(outfile_I, bmaj, bmin, bpa, '1', 'MWA', fname)
print 'Run aegean'

#cmd = 'mv %s STACKED' %(outfile_I)
#os.system(cmd)
