'''
Python program that reads in all images and then runs mirRegrid (miriad script) to regrid the images to match the first image (time-ordered). 
'''
import glob
import os

#######################################################################################
## MAIN PROGRAM       
#######################################################################################

##Read in the 2-min images and sort by name (time)
images = glob.glob('*MFS-image-pb.fits')
images.sort()

#Use first image as reference image
ref_image = images[0]


#For each image either copy (if reference image) or 
#run mirRegrid to regrid image to reference
for i in range(len(images)):
	if i > 0:
		img = images[i]
		path = os.getcwd()
		cmd = './mirRegrid %s %s' %(ref_image[:-5], img[:-5])
		os.system(cmd)
	else:
		new_img = ref_image[:-5]+'.regrid.fits'
		cmd = 'cp %s %s' %(ref_image, new_img)
		os.system(cmd)

##Clean up directory
os.system('rm -rf *.in')
os.system('rm -rf *regridded*')









