# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield

Significantly adapted from an earlier code by Farhad Allian.
"""

import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
import numpy as np
import scipy as sc
import scipy.ndimage
import gauss_fitting as gf
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")


class Full_map:
    def __init__(self, file_path, time_range=None, pixel_size=50.,
                 cadence=7.68):
        """
        Inputs:
            file_path = path to directory containing .fits files,
            time_range = slice(tinit, tfinal),
            pixel_size = size of pixels used in the observation (km),
            cadence = time between observation time frames (s).
        """
        self.file_path = os.path.abspath(file_path)
        hdul_list = []
        for i in range(238):
            hdul_list.append(fits.open(file_path + "/destretched_" + "%04d" % i + ".fits"))
        self.total_maps = hdul_list
        # set time range to whole time range if None
        if time_range is None:
            self.time_range = slice(0, len([name for name in os.listdir(self.file_path) if os.path.isfile(os.path.join(self.file_path, name))]))
        elif type(time_range) is slice:
            self.time_range = time_range
        else:
            raise ValueError("time_range must be of type slice")
        self.time_range = [range(0, 238)[self.time_range][0],
                           range(0, 238)[self.time_range][-1]]
        self.pixel_size = pixel_size
        self.cadence = cadence

    def crop(self, bottom_left, top_right):
        """
        Crop the map.

        Inputs:
            bottom_left = [x, y],
            top_right = [x, y].
        """
        for i, m in enumerate(self.total_maps):
            self.total_maps[i][0].data = m[0].data[bottom_left[0]:top_right[0],
                                                   bottom_left[1]:top_right[1]]

    def animate(self, slit_coords=None, interval=200, repeat_delay=1000,
                savefig=None):
        """
        Plot an animation of the data with optional overlayed slit.

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal],
            interval = time interval between frames in animation (ms),
            repeat_delay = time delay between repeats in animation (ms),
            savefig = None (not saved) or saved name string.
        """

        number_of_frames = self.time_range[1] - self.time_range[0]

        fig = plt.figure()

        im = plt.imshow(self.total_maps[self.time_range[0]][0].data,
                        aspect='auto', interpolation=None, animated=True,
                        origin='lower')

        def updatefig(i):
            image_data = self.total_maps[self.time_range[0] + i][0].data
            im.set_array(image_data)
            return im

        ani = animation.FuncAnimation(fig, updatefig, frames=number_of_frames,
                                      interval=interval,
                                      repeat_delay=repeat_delay)

        if slit_coords is not None:
            plt.plot(slit_coords[:2], slit_coords[2:], color='white')

        plt.show()

        if savefig is not None:
            ani.save(savefig)

    def distancetime(self, slit_coords, plot=False, savefig=None):
        """
        Distances in pixels

        Creates the distance-time plots along the slice coordiantes.

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal],
            savefig = None (not saved) or saved name string.
        """
        time_range_s = [self.time_range[0] * self.cadence,
                        self.time_range[1] * self.cadence]  # in s

        print("\nCalculating the intensity values along slit....")

        intensity1 = []

        # Number of points we can interpolate along the slit,
        # using pytagoras' theorem
        num = np.sqrt((slit_coords[1] - slit_coords[0])**2 + (slit_coords[3] - slit_coords[2])**2)
        print("num = " + str(num))

        x = np.linspace(slit_coords[0], slit_coords[1], num)
        y = np.linspace(slit_coords[2], slit_coords[3], num)

        for i, m in enumerate(self.total_maps):
            intensity1.append(sc.ndimage.map_coordinates(np.transpose(m[0].data),
                                                         np.vstack((x, y))))

        intensity1 = np.array(intensity1)

        if plot is True:
            plt.figure()
            plt.imshow(intensity1.T[:, self.time_range[0]:self.time_range[1]],
                       aspect='auto', interpolation=None, origin='lower',
                       extent=[time_range_s[0], time_range_s[1],
                               0, num*self.pixel_size])
            plt.xlabel('Time (s)')
            plt.ylabel('Distance (km)')
            plt.ylim([0, num*self.pixel_size])
            if savefig is not None:
                plt.savefig(savefig)

        return intensity1

    def find_boundaries(self, slit_coords, p0=[0.5, 45., 10., -1.1],
                        plot=False, savefig=None):
        """
        Find boundaries of the structure using gauss fitting. The edges are the
        half-maximum points on each side of the structure.

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal],
            p0 = initial [amplitude, mean, standard deviation, offset],
            savefig = None (not saved) or saved name string.
        """
        time_range_s = [self.time_range[0] * self.cadence,
                        self.time_range[1] * self.cadence]  # in s
        number_of_frames = self.time_range[1] - self.time_range[0]
        time_vals = np.linspace(self.time_range[0], self.time_range[1],
                                number_of_frames + 1)

        intensity1 = self.distancetime(slit_coords=slit_coords)

        boundary_t_vals = []
        boundary_x_vals_t = []
        boundary_x_vals_b = []

        for i, t in enumerate(time_vals):
            t = int(t)
            # Skip points which raise errors in gauss fitting.
            try:
                params = gf.gauss_fit(-intensity1[t], p0=p0,
                                      retrn="pos_half_max")
            except:
                pass
            boundary_t_vals.append(t * self.cadence)
            boundary_x_vals_b.append(params[0][0] * self.pixel_size)
            boundary_x_vals_t.append(params[0][1] * self.pixel_size)

        if plot is True:
            num = np.sqrt((slit_coords[1] - slit_coords[0])**2 + (slit_coords[3] - slit_coords[2])**2)
            plt.figure()
            plt.imshow(intensity1.T[:, self.time_range[0]:self.time_range[1]],
                       aspect='auto', interpolation=None, origin='lower',
                       extent=[time_range_s[0], time_range_s[1],
                               0, num*self.pixel_size])
            plt.xlabel('Time (s)')
            plt.ylabel('Distance (km)')
            plt.plot(boundary_t_vals, boundary_x_vals_b, 'wo')
            plt.plot(boundary_t_vals, boundary_x_vals_t, 'wo')
            plt.ylim([0, num*self.pixel_size])
            if savefig is not None:
                plt.savefig(savefig)

        return [boundary_x_vals_b, boundary_x_vals_t, boundary_t_vals]

    def intensity_slice(self, slit_coords, p0=[0.5, 45., 10., -1.1], gauss_fit=False,
                        savefig=None):
        """
        Create intensity slices for t in time_frames

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal],
            p0 = initial [amplitude, mean, standard deviation, offset],
            savefig = None (not saved) or saved name string.
        """

        for tf in range(self.time_range[0], self.time_range[1]):
            print("\nSlicing intensity at frame number " + str(tf))
            intensity1 = self.distancetime(slit_coords=slit_coords)
            intensity_slice = -intensity1[tf]
            s_vals = np.arange(len(intensity_slice))

            plt.figure()
            plt.plot(s_vals, intensity_slice, '-')

            if gauss_fit is True:
                # Plot the data with the best-fit model
                func = gf.gauss_fit(intensity_slice, p0=p0, retrn="func")
                boundaries = gf.gauss_fit(intensity_slice, p0=p0,
                                          retrn="pos_half_max")

                plt.plot(s_vals, func(s_vals), '-', color='red')
                plt.plot(boundaries[0], boundaries[1], 'ro')

                if savefig is not None:
                    plt.savefig(savefig)
