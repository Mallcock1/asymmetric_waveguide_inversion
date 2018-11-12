# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 10:32:22 2018

@author: Matt
"""

import sunpy.map
from astropy import units as u
from astropy.coordinates import SkyCoord

dir_path = 'D:/my_work/projects/Asymmetric_slab/Application_to_observation/Morton_2012_data/rosa_data'

all_images = sunpy.map.Map(dir_path + '/*fits', sequence=True)

#all_images.peek()


#total_maps = sunpy.map.Map((str(file_path)+'/*fits'),sequence=True)


## Define a region of interest
#length = 250
#x0 = -100
#y0 = -400

# Create a SunPy Map, and a second submap over the region of interest.
#smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
bottom_left = SkyCoord(-400*u.arcsec, -400*u.arcsec,
                       frame=all_images[0].coordinate_frame)
top_right = SkyCoord(358*u.arcsec, 294*u.arcsec,
                     frame=all_images[0].coordinate_frame)


submap1 = all_images[0].submap(bottom_left, top_right)