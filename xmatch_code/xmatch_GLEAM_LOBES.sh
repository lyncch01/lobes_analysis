##########################################################################################
# Match GLEAM with LOBES 093
##########################################################################################
START_TIME=$SECONDS
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
 	--table2=/home/corvus/Research/LOBES/catalogs/LOBES_chn093.fits \
 	--details2=NAMES,RA,ERR_RA,DEC,ERR_DEC,119,INT_FLUX,ERR_INT_FLUX,PA,ERR_PA,a,b,-,- \
 	--units2=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
 	--ra_lims2=0,360 \
 	--dec_lims2=-36,-18 \
 	--prefix2=lobes_093
STRING_1="Time to match LOBES 119 MHz: "
TIME_1=$SECONDS
ELAPSED_1=$(($TIME_1 - $START_TIME))
PRINT_1=$STRING_1$ELAPSED_1
echo $PRINT_1

##########################################################################################
# Match GLEAM with LOBES 121
##########################################################################################
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
 	--table2=/home/corvus/Research/LOBES/catalogs/LOBES_chn121.fits \
 	--details2=NAMES,RA,ERR_RA,DEC,ERR_DEC,154,INT_FLUX,ERR_INT_FLUX,PA,ERR_PA,a,b,-,- \
 	--units2=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
 	--ra_lims2=0,360 \
 	--dec_lims2=-36,-18 \
 	--prefix2=lobes_121
STRING_2="Time to match LOBES 154 MHz: "
TIME_2=$SECONDS
ELAPSED_2=$(($TIME_2 - $TIME_1))
PRINT_2=$STRING_2$ELAPSED_2
echo $PRINT_2

##########################################################################################
# Match GLEAM with LOBES 145
##########################################################################################
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
 	--table2=/home/corvus/Research/LOBES/catalogs/LOBES_chn121.fits \
 	--details2=NAMES,RA,ERR_RA,DEC,ERR_DEC,185,INT_FLUX,ERR_INT_FLUX,PA,ERR_PA,a,b,-,- \
 	--units2=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
 	--ra_lims2=0,360 \
 	--dec_lims2=-36,-18 \
 	--prefix2=lobes_145
STRING_3="Time to match LOBES 185 MHz: "
TIME_3=$SECONDS
ELAPSED_3=$(($TIME_3 - $TIME_2))
PRINT_3=$STRING_3$ELAPSED_3
echo $PRINT_3

##########################################################################################
# Match GLEAM with LOBES 169
##########################################################################################
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
 	--table2=/home/corvus/Research/LOBES/catalogs/LOBES_chn169.fits \
 	--details2=NAMES,RA,ERR_RA,DEC,ERR_DEC,216,INT_FLUX,ERR_INT_FLUX,PA,ERR_PA,a,b,-,- \
 	--units2=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
 	--ra_lims2=0,360 \
 	--dec_lims2=-36,-18 \
 	--prefix2=lobes_169
STRING_4="Time to match LOBES 216 MHz: "
TIME_4=$SECONDS
ELAPSED_4=$(($TIME_4 - $TIME_3))
PRINT_4=$STRING_4$ELAPSED_4
echo $PRINT_4

##########################################################################################
# Match GLEAM with VLSSr
##########################################################################################
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
	--table2=/home/corvus/Research/LOBES/catalogs/vlssr_names.fits \
	--details2=Name,RA_J2000,RA_err,DEC_J2000,DEC_err,74,Flux,Flux_err,PA,major,minor,-,- \
	--units2=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims2=0,140 \
	--dec_lims2=-20,40 \
	--prefix2=VLSSr
STRING_5="Time to match VLSSr: "
TIME_5=$SECONDS
ELAPSED_5=$(($TIME_5 - $TIME_4))
PRINT_5=$STRING_5$ELAPSED_5
echo $PRINT_5

##########################################################################################
# Match GLEAM with MRC
##########################################################################################
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
	--table2=/home/corvus/Research/LOBES/catalogs/vizier_mrc.fits \
	--details2=MRC,_RAJ2000,e_RA2000,_DEJ2000,e_DE2000,408,S408,e_S408,-,-,-,-,- \
	--units2=deg,sec,deg,arcsec,Jy,Jy,-,-,- \
	--ra_lims2=150,250 \
	--dec_lims2=-30,10 \
	--prefix2=MRC
STRING_6="Time to match MRC: "
TIME_6=$SECONDS
ELAPSED_6=$(($TIME_6 - $TIME_5))
PRINT_6=$STRING_6$ELAPSED_6
echo $PRINT_6

##########################################################################################
# Match GLEAM with SUMSS
##########################################################################################
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
	--table2=/home/corvus/Research/LOBES/catalogs/sumss_names.fits \
	--details2=name,_RAJ2000,e_RAJ2000,_DEJ2000,e_DEJ2000,843,St,e_St,PA,MajAxis,MinAxis,-,- \
	--units2=deg,arcsec,deg,arcsec,mJy,mJy,deg,arcsec,arcsec \
	--ra_lims2=0,100 \
	--dec_lims2=-60,-30 \
	--prefix2=SUMSS
STRING_7="Time to match SUMSS: "
TIME_7=$SECONDS
ELAPSED_7=$(($TIME_7 - $TIME_6))
PRINT_7=$STRING_7$ELAPSED_7
echo $PRINT_7

##########################################################################################
# Match GLEAM with NVSS
##########################################################################################
python $PUMA_DIR/scripts/cross_match.py --make_table --separation=140 \
 	--table1=/home/corvus/Research/LOBES/catalogs/GLEAM_EGC_v2.fits \
 	--details1=GLEAM,RAJ2000,err_RAJ2000,DEJ2000,err_DEJ2000,76,int_flux_076,err_int_flux_076,pa_076,a_076,b_076,-,-,84,int_flux_084,err_int_flux_084,92,int_flux_092,err_int_flux_092,99,int_flux_099,err_int_flux_099,107,int_flux_107,err_int_flux_107,115,int_flux_115,err_int_flux_115,122,int_flux_122,err_int_flux_122,130,int_flux_130,err_int_flux_130,143,int_flux_143,err_int_flux_143,151,int_flux_151,err_int_flux_151,158,int_flux_158,err_int_flux_158,166,int_flux_166,err_int_flux_166,174,int_flux_174,err_int_flux_174,181,int_flux_181,err_int_flux_181,189,int_flux_189,err_int_flux_189,197,int_197,err_int_flux_197,204,int_flux_204,err_int_flux_204,212,int_flux_212,err_int_flux_212,220,int_flux_220,err_int_flux_220,227,int_flux_227,err_int_flux_227 \
 	--units1=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
	--ra_lims1=135,175 \
	--dec_lims1=-30,0 \
	--prefix1=GLEAM \
	--table2=/home/corvus/Research/LOBES/catalogs/vizier_nvss.fits \
	--details2=NVSS,_RAJ2000,e_RAJ2000,_DEJ2000,e_DEJ2000,1400,S1_4,e_S1_4,PA,MajAxis,MinAxis,-,- \
	--units2=deg,sec,deg,arcsec,mJy,mJy,deg,arcsec,arcsec \
	--ra_lims2=100,200 \
	--dec_lims2=-40,0 \
	--prefix2=NVSS
STRING_8="Time to match NVSS: "
TIME_8=$SECONDS
ELAPSED_8=$(($TIME_8 - $TIME_7))
PRINT_8=$STRING_8$ELAPSED_8
echo $PRINT_8 

##########################################################################################
# Calculate Probabilities
##########################################################################################
python $PUMA_DIR/scripts/calculate_bayes.py --primary_cat=GLEAM --matched_cats=lobes_093,lobes_121,lobes_145,lobes_169,VLSSr,MRC,SUMSS,NVSS \
	--primary_freq=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227 --matched_freqs=119,154,185,216,74,408,843,1400 \
	--out_name=bayes_lobes.txt --resolution=00:02:20
STRING_9="Time to calculate posterior probabilities: "
TIME_9=$SECONDS
ELAPSED_9=$(($TIME_9 - $TIME_8))
PRINT_9=$STRING_9$ELAPSED_9
echo $PRINT_9

##########################################################################################
# Make Match Tables
##########################################################################################
python $PUMA_DIR/scripts/make_table.py --matched_cats=GLEAM,lobes_093,lobes_121,lobes_145,lobes_169,VLSSr,MRC,SUMSS,NVSS \
	--pref_cats=NVSS,SUMSS,lobes_169,lobes_145,lobes_121,lobes_093,GLEAM,MRC,VLSSR \
	--input_bayes=bayes_lobes.txt \
	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,119,154,185,216,74,408,843,1400 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:02:20 \
	--output_name=puma_lobes_split --verbose --split=00:02:00 \
	--big_cat --prefix=PUMA --format=fits
STRING_10="Time to apply criteria and create final matched table: "
TIME_10=$SECONDS
ELAPSED_10=$(($TIME_10 - $TIME_9))
PRINT_10=$STRING_10$ELAPSED_10
echo $PRINT_10

##########################################################################################
STRING_11="Total Time: "
TIME_11=$SECONDS
ELAPSED_11=$(($TIME_11 - $START_TIME))
PRINT_1=$STRING_11$ELAPSED_11
echo $PRINT_1
echo $PRINT_2
echo $PRINT_3
echo $PRINT_4
echo $PRINT_5
echo $PRINT_6
echo $PRINT_7
echo $PRINT_8
echo $PRINT_9
echo $PRINT_10
echo $PRINT_11

