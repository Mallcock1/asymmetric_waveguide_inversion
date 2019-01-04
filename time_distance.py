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
from astropy.convolution import convolve, Box1DKernel
# Nb: this overrides pylab's convolve
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
                        origin='lower', cmap="afmhot")

        def updatefig(i):
            image_data = self.total_maps[self.time_range[0] + i][0].data
            im.set_array(image_data)
            return im

        ani = animation.FuncAnimation(fig, updatefig, frames=number_of_frames,
                                      interval=interval,
                                      repeat_delay=repeat_delay)

        if slit_coords is not None:
            plt.plot(slit_coords[:2], slit_coords[2:], color='white')

            # testing moving average slits
            x0 = slit_coords[0]
            x1 = slit_coords[1]
            y0 = slit_coords[2]
            y1 = slit_coords[3]

            m_prime = (x0 - x1) / (y1 - y0)

            alpha = 5.  # / np.sqrt(1 + m_prime**2)
            for i in range(1, 3):
                plt.plot([x0 + i*alpha, x1 + i*alpha],
                         [y0 + i*alpha*m_prime, y1 + i*alpha*m_prime],
                         color='yellow')
                plt.plot([x0 - i*alpha, x1 - i*alpha],
                         [y0 - i*alpha*m_prime, y1 - i*alpha*m_prime],
                         color='green')

        plt.show()

        if savefig is not None:
            ani.save(savefig)

    def distancetime(self, slit_coords, moving_average=False, num_wtd_av=5,
                     wtd_av_distance=1., plot=False, savefig=None):
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

        intensity = []

        # Number of points we can interpolate along the slit,
        # using pytagoras' theorem
        num = np.sqrt((slit_coords[1] - slit_coords[0])**2
                      + (slit_coords[3] - slit_coords[2])**2)
        print("num = " + str(num))

        x = np.linspace(slit_coords[0], slit_coords[1], num)
        y = np.linspace(slit_coords[2], slit_coords[3], num)

#        x = np.arange(slit_coords[0], slit_coords[1])
#        y = np.arange(slit_coords[2], slit_coords[3])
        map_coords = scipy.ndimage.map_coordinates
        if moving_average is True:
            """
            At each time step, the sliced intensity is averaged over
            num_wtd_av number of pixels each side of the middle slit.

            This will increase the signal-to-noise ratio.

            num_wtd_av ~ 5 will give a good result.

            Higher will blur out large structures that might be interesting.

            High num_wtd_av is very computationally expensive.
            """
            x0 = slit_coords[0]
            x1 = slit_coords[1]
            y0 = slit_coords[2]
            y1 = slit_coords[3]

            m_prime = (x0 - x1) / (y1 - y0)

            alpha = wtd_av_distance  # / np.sqrt(1 + m_prime**2)

            for i, m in enumerate(self.total_maps):
                map_list_for_average = []
                # append middle slice
                map_list_for_average.append(map_coords(np.transpose(m[0].data),
                                            np.vstack((x, y))))
                for i in range(1, num_wtd_av + 1):
                    # append num_wtd_av number of adjacent slices
                    prev_slit = map_coords(np.transpose(m[0].data),
                                           np.vstack((x - i*alpha,
                                                      y - i*alpha*m_prime)))
                    next_slit = map_coords(np.transpose(m[0].data),
                                           np.vstack((x + i*alpha,
                                                      y - i*alpha*m_prime)))
                    map_list_for_average.append(prev_slit)
                    map_list_for_average.append(next_slit)
                # take average of these slices
                m_average = list(np.mean(map_list_for_average, axis=0))
                intensity.append(m_average)

        else:
            for i, m in enumerate(self.total_maps):
                intensity.append(map_coords(np.transpose(m[0].data),
                                             np.vstack((x, y))))

        intensity = np.array(intensity)

        if plot is True:
            plt.figure()
            plt.imshow(intensity.T[:, self.time_range[0]:self.time_range[1]],
                       aspect='auto', interpolation=None, origin='lower',
                       extent=[time_range_s[0], time_range_s[1],
                               0, num*self.pixel_size], cmap="afmhot")
            plt.xlabel('Time (s)')
            plt.ylabel('Distance (km)')
            plt.ylim([0, num*self.pixel_size])
            if savefig is not None:
                plt.savefig(savefig)

        return intensity

    def find_boundaries(self, slit_coords, moving_average=False,
                        wtd_av_distance=1., num_wtd_av=5,
                        p0=[0.5, 45., 10., -1.1], stabilise=False,
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

        intensity = self.distancetime(slit_coords=slit_coords,
                                      moving_average=moving_average,
                                      wtd_av_distance=wtd_av_distance,
                                      num_wtd_av=num_wtd_av)

        boundary_t_vals = []
        boundary_x_vals_t = []
        boundary_x_vals_b = []

        FWHM_factor = np.sqrt(2*np.log(2))
        success = False
        for i, t in enumerate(time_vals):
            t = int(t)
            if stabilise is False:
                data_to_fit = -intensity[t]
                if p0[1] is None:
                    # set initial guess at gaussian mean to be argmax of 
                    # intensity
                    p0[1] = np.argmax(data_to_fit)
                # Skip points which raise errors in gauss fitting.
                try:
                    params = gf.gauss_fit(data_to_fit, p0=p0, retrn="params")
                except RuntimeError:
                    pass
    
                # bottom and top x_vals
                bot = (params[1] - np.sqrt(2*np.log(2))*params[2])*self.pixel_size
                top = (params[1] + np.sqrt(2*np.log(2))*params[2])*self.pixel_size
    
                # Just append the first one
                if i == 0:
                    boundary_t_vals.append(t * self.cadence)
                    boundary_x_vals_b.append(bot)
                    boundary_x_vals_t.append(top)
                else:
                    # if big jump in width, just skip and use prev params for next
                    if top - bot < 2*(boundary_x_vals_t[-1] - boundary_x_vals_b[-1]):
                        boundary_t_vals.append(t * self.cadence)
                        boundary_x_vals_b.append(bot)
                        boundary_x_vals_t.append(top)
    
                        p0 = params
            else:
                ###############################
                # Crop gaussian range if jump
                ###############################
                data_to_fit = -intensity[t]
                if p0[1] is None:
                    # set initial guess at gaussian mean to be argmax of 
                    # intensity
                    p0[1] = np.argmax(data_to_fit)
                # Skip points which raise errors in gauss fitting.
                try:
                    params = gf.gauss_fit(data_to_fit, p0=p0, retrn="params")
                except RuntimeError:
                    pass
    
                # bottom and top x_vals
                bot = (params[1] - np.sqrt(2*np.log(2))*params[2])*self.pixel_size
                top = (params[1] + np.sqrt(2*np.log(2))*params[2])*self.pixel_size
    
                # Just append the first one
                if i == 0:
                    boundary_t_vals.append(t * self.cadence)
                    boundary_x_vals_b.append(bot)
                    boundary_x_vals_t.append(top)
                else:
                    # if big jump in width, just skip and use prev params for next
                    if top - bot < 2*(boundary_x_vals_t[-1] - boundary_x_vals_b[-1]):
                        boundary_t_vals.append(t * self.cadence)
                        boundary_x_vals_b.append(bot)
                        boundary_x_vals_t.append(top)
    
                        p0 = params
                    else:
                        # crop distance range to just around previous boundaries
                        bot_boundary_prev = int(np.round(p0[1] - FWHM_factor*p0[2]))
                        top_boundary_prev = int(np.round(p0[1] + FWHM_factor*p0[2]))
                        
                        bot_of_data = max(bot_boundary_prev - stabilise, 0)
                        top_of_data = min(top_boundary_prev + stabilise,
                                          len(-intensity[t]))
    
                        data_to_fit = -intensity[t][bot_of_data:top_of_data]
    
                        p0[1] = p0[1] - bot_of_data
                        try:
                            params = gf.gauss_fit(data_to_fit, p0=p0,
                                                  retrn="params")
    
                            p0_new = params
                            p0_new[1] = p0_new[1] + bot_of_data
    
                            bot = (p0_new[1] - FWHM_factor*p0_new[2])*self.pixel_size
                            top = (p0_new[1] + FWHM_factor*p0_new[2])*self.pixel_size
    
                            # if big jump in width, just skip and use prev params for next
                            if top - bot < 2.5*(boundary_x_vals_t[-1] - boundary_x_vals_b[-1]):
                                boundary_t_vals.append(t * self.cadence)
                                boundary_x_vals_b.append(bot)
                                boundary_x_vals_t.append(top)
    
                                p0 = p0_new
                        except RuntimeError:
                            pass
                
                
                
                
                
                
                
                
#                if success is False:
#                    data_to_fit = -intensity[t]
#                    # Skip points which raise errors in gauss fitting.
#                    try:
#                        if p0[1] is None:
#                            # set initial guess at gaussian mean to be argmax of 
#                            # intensity
#                            p0[1] = np.argmax(data_to_fit)
#                        params = gf.gauss_fit(data_to_fit, p0=p0,
#                                              retrn="params")
#                        success = True
#    
#                        # bottom and top x_vals
#                        bot = (params[1] - FWHM_factor*params[2])*self.pixel_size
#                        top = (params[1] + FWHM_factor*params[2])*self.pixel_size
#    
#                        boundary_t_vals.append(t * self.cadence)
#                        boundary_x_vals_b.append(bot)
#                        boundary_x_vals_t.append(top)
#                        
#                        p0 = params
#                    except RuntimeError:
#                        pass
#                else:
                    ###############################
                    # Crop gaussian range if jump
                    ###############################
#                    # crop distance range to just around previous boundaries
#                    bot_boundary_prev = int(np.round(p0[1] - FWHM_factor*p0[2]))
#                    top_boundary_prev = int(np.round(p0[1] + FWHM_factor*p0[2]))
#                    
#                    bot_of_data = max(bot_boundary_prev - stabilise, 0)
#                    top_of_data = min(top_boundary_prev + stabilise,
#                                      len(-intensity[t]))
#
#                    data_to_fit = -intensity[t][bot_of_data:top_of_data]
#
#                    p0[1] = p0[1] - bot_of_data
#                    try:
#                        params = gf.gauss_fit(data_to_fit, p0=p0,
#                                              retrn="params")
#
#                        p0_new = params
#                        p0_new[1] = p0_new[1] + bot_of_data
#
#                        bot = (p0_new[1] - FWHM_factor*p0_new[2])*self.pixel_size
#                        top = (p0_new[1] + FWHM_factor*p0_new[2])*self.pixel_size
#
#                        # if big jump in width, just skip and use prev params for next
#                        if top - bot < 2.5*(boundary_x_vals_t[-1] - boundary_x_vals_b[-1]):
#                            boundary_t_vals.append(t * self.cadence)
#                            boundary_x_vals_b.append(bot)
#                            boundary_x_vals_t.append(top)
#
#                            p0 = p0_new
#                    except RuntimeError:
#                        pass

        if plot is True:
            num = np.sqrt((slit_coords[1] - slit_coords[0])**2
                          + (slit_coords[3] - slit_coords[2])**2)
            plt.figure()
            plt.imshow(intensity.T[:, self.time_range[0]:self.time_range[1]],
                       aspect='auto', interpolation=None, origin='lower',
                       extent=[time_range_s[0], time_range_s[1],
                               0, num*self.pixel_size], cmap="afmhot")
            plt.xlabel('Time (s)')
            plt.ylabel('Distance (km)')
            plt.plot(boundary_t_vals, boundary_x_vals_b, 'wo')
            plt.plot(boundary_t_vals, boundary_x_vals_t, 'wo')
            plt.ylim([0, num*self.pixel_size])
            if savefig is not None:
                plt.savefig(savefig)

        return [boundary_x_vals_b, boundary_x_vals_t, boundary_t_vals]

    def intensity_slice(self, slit_coords, time_slice,
                        p0=[0.5, 45., 10., -1.1], gauss_fit=False,
                        savefig=None):
        """
        Create intensity slices for t in time_frames

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal],
            p0 = initial [amplitude, mean, standard deviation, offset],
            savefig = None (not saved) or saved name string.
        """

        for ts in time_slice:
            print("\nSlicing intensity at frame number " + str(ts))
            intensity = self.distancetime(slit_coords=slit_coords)
            intensity_slice = -intensity[ts]
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

    def find_multi_slit_boundaries(self, slit_coords, num_slits=1,
                                   slit_distance=1., moving_average=False,
                                   wtd_av_distance=1., num_wtd_av=5,
                                   p0=[0.5, 45., 10., -1.1], stabilise=False,
                                   plot=False, savefig=None):
        """
#        Find boundaries of the structure using gauss fitting. The edges are the
#        half-maximum points on each side of the structure.

        Inputs:
            slit_coords = [xinit, xfinal, yinit, yfinal],
            p0 = initial [amplitude, mean, standard deviation, offset],
            savefig = None (not saved) or saved name string.
        """

        if num_slits % 2 != 1:
            raise ValueError("num_slits must be an odd number.")

        multi_boundaries = []
        middle_slit = self.find_boundaries(slit_coords=slit_coords,
                                           moving_average=moving_average,
                                           wtd_av_distance=wtd_av_distance,
                                           num_wtd_av=num_wtd_av, p0=p0,
                                           stabilise=stabilise)
        multi_boundaries.append(middle_slit)

        x0 = slit_coords[0]
        x1 = slit_coords[1]
        y0 = slit_coords[2]
        y1 = slit_coords[3]
        m_prime = (x0 - x1) / (y1 - y0)
        alpha = slit_distance  # / np.sqrt(1 + m_prime**2)
        for i in range(1, int((num_slits + 1) / 2)):
            x0_r = x0 + i*alpha
            x1_r = x1 + i*alpha
            y0_r = y0 + i*alpha*m_prime
            y1_r = y1 + i*alpha*m_prime
            slit_coords_r = [x0_r, x1_r, y0_r, y1_r]

            x0_l = x0 - i*alpha
            x1_l = x1 - i*alpha
            y0_l = y0 - i*alpha*m_prime
            y1_l = y1 - i*alpha*m_prime
            slit_coords_l = [x0_l, x1_l, y0_l, y1_l]

            slit_r = self.find_boundaries(slit_coords=slit_coords_r,
                                          moving_average=moving_average,
                                          wtd_av_distance=wtd_av_distance,
                                          num_wtd_av=num_wtd_av, p0=p0,
                                          stabilise=stabilise)

            slit_l = self.find_boundaries(slit_coords=slit_coords_l,
                                          moving_average=moving_average,
                                          wtd_av_distance=wtd_av_distance,
                                          num_wtd_av=num_wtd_av, p0=p0,
                                          stabilise=stabilise)
            multi_boundaries.append(slit_r)
            multi_boundaries.insert(0, slit_l)

        if plot is True:
            num = np.sqrt((slit_coords[1] - slit_coords[0])**2
                          + (slit_coords[3] - slit_coords[2])**2)
            plt.figure()
            plt.xlabel('Time (s)')
            plt.ylabel('Distance (km)')
            plt.ylim([0, num*self.pixel_size])
            for i, boundary in enumerate(multi_boundaries):
                boundary_shift_b = np.array(boundary[0]) + i*slit_distance*self.pixel_size
                boundary_shift_t = np.array(boundary[1]) + i*slit_distance*self.pixel_size
                plt.plot(boundary[2], boundary_shift_b, 'o')
                plt.plot(boundary[2], boundary_shift_t, 'o')
                plt.plot(boundary[2], convolve(boundary_shift_b, Box1DKernel(3)))
                plt.plot(boundary[2], convolve(boundary_shift_t, Box1DKernel(3)))
            if savefig is not None:
                plt.savefig(savefig)

        return multi_boundaries

    def find_multi_slit_widths(self, slit_coords, num_slits=1, slit_distance=1.,
                               moving_average=False, wtd_av_distance=1.,
                               num_wtd_av=5, p0=[0.5, 45., 10., -1.1], stabilise=False,
                               plot=False, savefig=None):
        multi_boundaries = self.find_multi_slit_boundaries(slit_coords=slit_coords,
                                                           num_slits=num_slits,
                                                           slit_distance=slit_distance,
                                                           moving_average=moving_average,
                                                           wtd_av_distance=wtd_av_distance,
                                                           num_wtd_av=num_wtd_av,
                                                           p0=p0, stabilise=stabilise)
        widths = []
        for i, mb in enumerate(multi_boundaries):
            widths.append([np.array(mb[1]) - np.array(mb[0]),
                           np.array(mb[2])])

        if plot is True:
#            num = np.sqrt((slit_coords[1] - slit_coords[0])**2
#                          + (slit_coords[3] - slit_coords[2])**2)
            plt.figure()

            plt.xlabel('Time (s)')
            plt.ylabel('Distance (km)')
#            plt.ylim([0, num*self.pixel_size])
            color_list = [(0,0,0), (0.25,0,0), (0.5,0,0), (0.75,0,0), (1,0,0)]
            for i, width in enumerate(widths):
                width_shift = width[0] + i*slit_distance*self.pixel_size
                plt.errorbar(width[1], width_shift, yerr=self.pixel_size,
                             fmt='none', color=color_list[i])
                plt.plot(width[1], convolve(width_shift, Box1DKernel(3)),
                         color=color_list[i])
            if savefig is not None:
                plt.savefig(savefig)
        return widths
