# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 08:58:12 2018

@author: Matt
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits


########################################
#Pickling the data
dir_path = 'D:/my_work/projects/Asymmetric_slab/Application_to_observation/Morton_2012_data/rosa_data/'
file_path = dir_path + 'destretched_0000.fits'

#display the structure of the file:
fits.info(file_path)

image_time_data = fits.getdata(file_path)[:720, :780]

# read the image data
for i in range(1,238):
    file_path = dir_path + 'destretched_0' + "%03d" % i + '.fits'
    image_time_data = np.dstack((image_time_data, fits.getdata(file_path)[:720, :780]))
    
import pickle
pickling_on = open("image_time_data.pickle","wb")
pickle.dump(image_time_data, pickling_on)
pickling_on.close()
##########################################

#unpickle the data
import pickle
pickle_off = open("image_time_data.pickle","rb")
image_time_data = pickle.load(pickle_off)

#########################################
#
##Display the image data:
#plt.figure()
#plt.imshow(image_time_data[:,:,23])
#plt.colorbar()


########################################
#animate through time
fig = plt.figure()

#animation
im = plt.imshow(image_time_data[:,:,0], animated=True)

#i = 1
def updatefig(i):
    image_data = image_time_data[:,:,i]
    im.set_array(image_data)
    return im

ani = animation.FuncAnimation(fig, updatefig, frames=len(image_time_data[0,0,:]), interval=1)

########################################

