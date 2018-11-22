# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield

Adapted from an earlier code by Farhad Allian
"""

import matplotlib.pyplot as plt
import sunpy
import os
from matplotlib import colors
from astropy.coordinates import SkyCoord
import matplotlib.animation as animation
import astropy.units as u
from skimage import exposure
import numpy as np
import scipy as sc
import warnings
warnings.filterwarnings("ignore")


class Full_map:
    def __init__(self, file_path):
        self.file_path = os.path.abspath(file_path)

    def listdir_fullpath(d):
        return [os.path.join(d, f) for f in os.listdir(d)]

    total_files = listdir_fullpath(location)
    total_files = sorted(total_files)

    if os.path.exists(file_path):
        print("\nDirectory exists.. ")
        #checking
        counter = 0
        for file in total_files:
            if file.endswith("s"):
                counter += 1

    if len(total_files) == counter:
        print("\nTotal number of FITS files in this directory is %s. \n" % len(total_files))


    def createmapcube(self):
        """
        This function creates the data/map cube.
    
        Inputs:
            getmaps- whether you want to create the mapcube in the given directory.
            plot- creates the standard plot of the first image.
            basic_plot- same as above, but without axes and labels.
            draw_grid- plots the solar longitude and latitude and the limb.
            draw_limb- plots only the limb.
            submap - whether you want to look at a specific region of interest.
            vmin,vmax (optional) - intensity of images.
    
        Example call:
            >> [in]: from dtp import createmapcube
            >> [in]: createmapcube(getmaps=True)
        """
        print("\nCreating datacube... ".format(len(total_files)))
        return sunpy.map.Map(str(self.file_path) + '/*ts', sequence=True)
    
        if len(total_maps) == len(total_files):
            print("\nMapcube created.")
        else:
            raise Exception("The number of files does not equal the number of "
                            "maps. Please create mapcube again.")


    total_maps = createmapcube(file_path)

    def distancetime(self, xfinal, xinitial, yfinal, yinitial):
        """
        Distances in pixels

        This function creates the distance-time plots along the slice coordiantes.

        Inputs:
        getintens- whether to calculate the intensity values along the slice.
        dt- whether you want to plot only distance time plot.
        peaks, threhsold, box_size- finds the maximum intensity value on distance 
                                    time plots. (For time series analysis).
        intensityslice - I = I(t,s(x,y)).
        """

        print("\nCalculating the intensity values along slit....")

        intensity1 = []

        # Number of points we can interpolate along the slit.
        num = np.sqrt((xfinal - xinitial)**2 + (yfinal - yinitial)**2)
        print("num = " + str(num))
        
        global x, y
        x = np.linspace(xfinal, xinitial, num)
        y = np.linspace(yfinal, yinitial, num)
    
        for i, m in enumerate(total_maps):
            intensity1.append(sc.ndimage.map_coordinates(np.transpose(m.data),
                                                         np.vstack((x, y))))
        cad = 7.68
        global time
        time = np.linspace(0, len(intensity1)*cad, len(intensity1))
    
        global space
        space = np.linspace(0, num*50., len(intensity1[0]))  # multiply by 50 to convert into km
    
        return np.array(intensity1)

#################
#        OPP EDIT FROM HERE

    time_window = time_window_frames * 7.68  # in s
    number_of_frames = time_window_frames[1] - time_window_frames[0]
    
    intensity1 = distancetime(xfinal=slit_coords_x[0], xinitial=slit_coords_x[1],
                              yfinal=slit_coords_y[0], yinitial=slit_coords_y[1])

print("\nSlicing intensity....")
intensity_slice = -intensity1[181]
plt.figure()
plt.plot(intensity_slice)



time_vals = np.linspace(time_window_frames[0], time_window_frames[1],
                        number_of_frames + 1)
# s is the coordinate along the slit
s_vals = np.arange(len(intensity1))  # num - 1

#fit a Gaussian curve
from scipy.optimize import curve_fit
def gauss_func(x, A, x0, sigma, offset): #a = amplitude, b = angular frequency
    return np.abs(A) * np.exp(-(x - x0)**2 / (2*sigma**2)) + offset

boundary_t_vals = []
boundary_x_vals_t = []
boundary_x_vals_b = []
#
for t in np.linspace(time_window_frames[0], time_window_frames[1] - 1, number_of_frames):
    t = int(t)
#    print('t = ' + str(t))
    params, params_covariance = curve_fit(gauss_func, s_vals, -intensity1[t], p0=[0.5, 45., 10., -1.1])
#    print('params =' + str(params))
    boundary_t_vals.append(t * 7.68)
    boundary_x_vals_t.append((params[1] + np.sqrt(2*np.log(2))*params[2]) * 50)
    boundary_x_vals_b.append((params[1] - np.sqrt(2*np.log(2))*params[2]) * 50)


    def td_gauss_plot(self.savefig=None):
        intensity1 = self.distancetime(xfinal=slit_coords_x[0],
                                       xinitial=slit_coords_x[1],
                                       yfinal=slit_coords_y[0],
                                       yinitial=slit_coords_y[1])
        plt.figure()
        plt.imshow(intensity1.T[:, time_window_frames[0]:time_window_frames[1]],
                   aspect='auto', interpolation=None, origin='lower',
                   extent=[time_window[0], time_window[1],
                           space.min(), space.max()])  # 1100, 1500 #To change x/y lims: add kwarg extent=[time_window[0], time_window[1], space.min(), space.max()]
        plt.xlabel('Time (s)')
        plt.ylabel('Distance (km)')
        plt.plot(boundary_t_vals, boundary_x_vals_b, 'wo')
        plt.plot(boundary_t_vals, boundary_x_vals_t, 'wo')
        plt.ylim([space.min(), space.max()])
        if savefig is not None:
            plt.savefig(savefig)








## Plot the data with the best-fit model
#params, params_covariance = curve_fit(gauss_func, x_vals, intensity_slice, p0=[0.5, 45., 10., -1.1])
#plt.figure()
#plt.plot(x_vals, intensity_slice, '-')
#plt.plot(x_vals, gauss_func(x_vals, params[0], params[1], params[2], params[3]), '-', color='red')
#plt.plot([params[1] - np.sqrt(2*np.log(2))*params[2], params[1] + np.sqrt(2*np.log(2))*params[2]],
#         [params[3] + 0.5*params[0], params[3] + 0.5*params[0]], 'ro')















######################################

#plt.figure()
#plt.imshow(total_maps[0].data, aspect='auto', interpolation=None, origin='lower')
#plt.plot(slit_coords_x, slit_coords_y, color='white')




##################################
#Animate through time
import matplotlib.animation as animation

fig = plt.figure()

def updatefig(i):
    image_data = total_maps[i].data
    im.set_array(image_data)
    return im

im = plt.imshow(total_maps[time_window_frames[0]].data, aspect='auto', 
                interpolation=None, animated=True, origin='lower')

ani = animation.FuncAnimation(fig, updatefig, frames=number_of_frames, interval=200, repeat_delay = 1000) #interval = delay between each fram in ms
plt.plot(slit_coords_x, slit_coords_y, color='white')


