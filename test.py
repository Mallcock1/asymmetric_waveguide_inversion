# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield
"""

import time_distance as td
import astropy.units as u


file_path = 'D:/my_work/projects/Asymmetric_slab/Application_to_observation/Morton_2012_data/rosa_data'

#morton12 = td.Full_map(file_path)

#morton12.distancetime(slit_coords=[100, 200, 100, 200], time_range=[0,238],
#                      plot=True, savefig=None)
#morton12.animate(slit_coords=[100, 200, 100, 200])
#morton12.find_boundaries(slit_coords=[100, 200, 100, 200], time_range=[0,238],
#                        p0=[0.5, 45., 10., -1.1], plot=True, savefig=None)

from astropy.io import fits

hdul_list = []

for i in range(238):
    hdul_list[i] = fits.open(file_path + "/destretched_" + "%04d" % i + ".fits")
    hdul_list[i][0].data = hdul_list[i][0].data[:720, :780]

import matplotlib.pyplot as plt

plt.imshow(hdul_list[i][0].data, aspect='auto', interpolation=None, origin='lower')

fig = plt.figure()

im = plt.imshow(self.total_maps[time_range[0]].data, aspect='auto',
                interpolation=None, animated=True, origin='lower')

def updatefig(i):
    image_data = self.total_maps[time_range[0] + i].data
    im.set_array(image_data)
    return im

ani = animation.FuncAnimation(fig, updatefig, frames=number_of_frames,
                              interval=200, repeat_delay=1000)  # interval = delay between each fram in ms

if slit_coords is not None:
    plt.plot(slit_coords[:2], slit_coords[2:], color='white')
    print(slit_coords)
plt.show(ani)

#
#
#morton12.crop([0, 100], [0, 500])
#
#morton12.cropped_maps.peek()

#morton12.total_maps.peek()
#
#morton12.animate()

#morton12.total_maps[0].pixel_to_world(100*u.pix, 0*u.pix, origin=1000)

#morton12.total_maps[0].top_right_coord

#morton12.total_maps[0].peek()
