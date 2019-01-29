# lobes_analysis
Code associated with analysis of LOBES data

Processing of LOBES data includes the following steps:

(1) Creation of median-stacked image from individual 2-min snapshot images.

For each LOBES field and frequency first use mirRegrid and regrid_images.py to regrid each individual 2-min image to a reference obsID (mirRegrid is run using miriad so be sure to have it installed). 

(2) Regridded images are then copied to Magnus to be combined to create the median stacked image. This is done using process_images.py and submitted to the queue using qprocessImages.sh. The python code process_images.py also runs BANE to create background/rms images and then does source finding using aegean (output is image_name_comp.fits). 

The median stacked images plus its background/rms images and the source catalog (comp.fits) are copied to my personal computer for further analysis. 

(3) Once source catalogs are created for each of the fields (either those around EoR1 or EoR0) for a single frequency, final_catalog_single.py is used to create a single catalog of sources that combines the sources from all the fields while removing sources repeated in the individual field catalogs (ensures sources are not repeated in the final catalog). This script only selects sources from the individual catalogs that are within the field-of-view for each frequency: 

119 MHz ==> 32 

154 MHz ==> 24 deg

185 MHz ==> 20 deg

216 MHz ==> 17 deg

(4) Using full source catalogs for each frequency, these are then cross-matched to create a single catalog with all spectral information for each source found in the survey; this is done using spec_cross.py. Sources only found at one frequency are recorded in this catalog. 

(5) PUMA is then run cross-matching the final LOBES spectral catalog to GLEAM using cross_match_GLEAM_LOBES.sh

