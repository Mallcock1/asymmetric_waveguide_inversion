# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield
"""

import time_distance as td

file_path = 'D:/my_work/projects/Asymmetric_slab/Application_to_observation/Morton_2012_data/rosa_data'

morton12 = td.Full_map(file_path)
morton12.crop([0,0], [720, 780])
#morton12.animate(slit_coords=[425, 429, 646, 546], interval=100, savefig="/plots/animation.mp4")
morton12.distancetime(slit_coords=[591, 580, 265, 371], time_range=[0,71],
                      plot=True, savefig="plots/test.png")
#morton12.intensity_slice(slit_coords=[591, 580, 265, 371], time_frames=[60],
#                         gauss_fit=True, savefig=None)
#boundaries = morton12.find_boundaries(slit_coords=[591, 580, 265, 371],
#                         time_range=[0,71], p0=[0.5, 45., 10., -1.1],
#                         plot=True, savefig=None)

#slit_coords_x = [591, 580]
#slit_coords_y = [265, 371]

#morton12.crop([0, 100], [0, 500])
#
#morton12.cropped_maps.peek()

#morton12.total_maps.peek()
#
#morton12.animate()

#morton12.total_maps[0].pixel_to_world(100*u.pix, 0*u.pix, origin=1000)

#morton12.total_maps[0].top_right_coord

#morton12.total_maps[0].peek()
