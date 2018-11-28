# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield

Adapted from an earlier code by Farhad Allian


"""

import matplotlib.pyplot as plt
import sunpy
import os
import matplotlib.animation as animation
import numpy as np
import scipy as sc
import gauss_fitting as gf
import warnings
warnings.filterwarnings("ignore")


class Full_map:
    def __init__(self, file_path):
        self.file_path = os.path.abspath(file_path)
        self.total_maps = sunpy.map.Map(str(self.file_path) + '/*.fits',
                                        sequence=True)

    def distancetime(self, slit_coords):
        """
        Distances in pixels

        Creates the distance-time plots along the slice coordiantes.

        Inputs:
        slit_coords = [xinit, xfinal, yinit, yfinal]
        """

        print("\nCalculating the intensity values along slit....")

        intensity1 = []

        # Number of points we can interpolate along the slit.
        num = np.sqrt((slit_coords[1] - slit_coords[0])**2 + (slit_coords[3] - slit_coords[2])**2)
        print("num = " + str(num))

        global x, y
        x = np.linspace(slit_coords[0], slit_coords[1], num)
        y = np.linspace(slit_coords[2], slit_coords[3], num)
    
        for i, m in enumerate(self.total_maps):
            intensity1.append(sc.ndimage.map_coordinates(np.transpose(m.data),
                                                         np.vstack((x, y))))
        cad = 7.68
        global time
        time = np.linspace(0, len(intensity1)*cad, len(intensity1))
    
        global space
        space = np.linspace(0, num*50., len(intensity1[0]))  # multiply by 50 to convert into km

        return np.array(intensity1)

    def find_boundaries(self, slit_coords, time_range,
                        p0=[0.5, 45., 10., -1.1], plot=False, savefig=None):
        """
        Find boundaries of the structure using gauss fitting. The edges are the
        half-maximum points on each side of the structure.

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal]
            time_range = [tinit, tfinal]
        """
        time_range_s = time_range * 7.68  # in s
        number_of_frames = time_range[1] - time_range[0]
        time_vals = np.linspace(time_range[0], time_range[1],
                                number_of_frames + 1)

        intensity1 = self.distancetime(slit_coords=slit_coords)

        boundary_t_vals = []
        boundary_x_vals_t = []
        boundary_x_vals_b = []

        for i, t in enumerate(time_vals):  # np.linspace(time_frames[0], time_frames[1] - 1, number_of_frames):
            t = int(t)
            params = gf.gauss_fitting(-intensity1[t], p0=p0, retrn="params")
            boundary_t_vals.append(params[2] * 7.68)
            boundary_x_vals_t.append(params[1] * 50)
            boundary_x_vals_b.append(params[2] * 50)

        return [boundary_x_vals_b, boundary_x_vals_t, boundary_t_vals]

        if plot is True:
            plt.figure()
            plt.imshow(intensity1.T[:, time_range[0]:time_range[1]],
                       aspect='auto', interpolation=None, origin='lower',
                       extent=[time_range_s[0], time_range_s[1],
                               space.min(), space.max()])  # 1100, 1500 #To change x/y lims: add kwarg extent=[time_window[0], time_window[1], space.min(), space.max()]
            plt.xlabel('Time (s)')
            plt.ylabel('Distance (km)')
            plt.plot(boundary_t_vals, boundary_x_vals_b, 'wo')
            plt.plot(boundary_t_vals, boundary_x_vals_t, 'wo')
            plt.ylim([space.min(), space.max()])
            if savefig is not None:
                plt.savefig(savefig)

    def intensity_slice(self, slit_coords, time_frames,
                        p0=[0.5, 45., 10., -1.1], gauss_fit=False):
        """
        Create intensity slices for t in time_frames

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal]
            time_frames = [t1, t2, ...]
        """

        for tf in time_frames:
            print("\nSlicing intensity at frame number " + str(tf))
            intensity1 = self.distancetime(slit_coords=slit_coords)
            intensity_slice = -intensity1[tf]
            s_vals = np.arange(len(intensity_slice))

            plt.figure()
            plt.plot(s_vals, intensity_slice, '-')

            if gauss_fit is True:
                # Plot the data with the best-fit model
                func = gf.gauss_fitting(intensity_slice, p0=p0, retrn="func")
                boundaries = gf.gauss_fitting(intensity_slice, p0=p0,
                                              retrn="pos_half_max")

                plt.plot(s_vals, func(s_vals), '-', color='red')
                plt.plot(boundaries[0], boundaries[1], 'ro')

    def animate(self, time_range=[0, 237], slit_coords=None):
        """
        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal]
            time_range = [tinit, tfinal]
        """
        fig = plt.figure()
        # ims is a list of lists, each row is a list of artists to draw in the
        # current frame; here we are just animating one artist, the image, in
        # each frame
        ims = []
        for i in range(time_range[0], time_range[1]):
            im = plt.imshow(self.total_maps[i].data, aspect='auto',
                            interpolation=None, animated=True, origin='lower')
            ims.append([im])
        
        ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                        repeat_delay=1000)
        
        # ani.save('dynamic_images.mp4')
        
        plt.show()
        
        
        

#        number_of_frames = time_range[1] - time_range[0]
#
#        fig = plt.figure()
#
#        im = plt.imshow(self.total_maps[time_range[0]].data, aspect='auto',
#                        interpolation=None, animated=True, origin='lower')
#
#        def updatefig(i):
#            image_data = self.total_maps[time_range[0] + i].data
#            im.set_array(image_data)
#            return im
#
#        ani = animation.FuncAnimation(fig, updatefig, frames=number_of_frames,
#                                      interval=200, repeat_delay=1000)  # interval = delay between each fram in ms
#        plt.show()
#
#        if slit_coords is not None:
#            plt.plot(slit_coords[:2], slit_coords[2:], color='white')
        
######################################

#plt.figure()
#plt.imshow(total_maps[0].data, aspect='auto', interpolation=None, origin='lower')
#plt.plot(slit_coords_x, slit_coords_y, color='white')