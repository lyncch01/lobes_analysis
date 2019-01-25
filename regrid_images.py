'''
Program runs within CASA environment and does the following: (1) create sub-images centred 
on targeted sources (coord read from cvs file); (2) fit 2min sub-images; (3) median stack 
images on specified time-scales for all sources; (4) fit median-stacked images; 
(5) create summary of all >4 sigma detections (median stack + snap shot)
'''

import numpy as np
import glob
from astropy.io import fits
from astropy import units as au
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import os
import gc
import csv




#######################################################################################
#								 MAIN PROGRAM       
#######################################################################################

#os.system("mkdir Images")
images = glob.glob('*MFS-image-pb.fits')
images.sort()
ref_image = images[0]
print images

for i in range(len(images)):
	if i > 0:
		img = images[i]
		path = os.getcwd()
		cmd = './mirRegrid %s %s' %(ref_image[:-5], img[:-5])
		print cmd
		os.system(cmd)
		#out_name = regridIm(images[i], ref_image)
		#print out_name
	else:
		new_img = ref_image[:-5]+'.regrid.fits'
		cmd = 'cp %s %s' %(ref_image, new_img)
		os.system(cmd)









