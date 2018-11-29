# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield
"""

import time_distance as td

file_path = 'D:/my_work/projects/Asymmetric_slab/Application_to_observation/Morton_2012_data/rosa_data'

morton12 = td.Full_map(file_path)

morton12.distancetime(slit_coords=[100, 200, 100, 200], time_range=[0,238],
                      plot=True, savefig=None)
morton12.animate(slit_coords=[100, 200, 100, 200])
#morton12.find_boundaries(slit_coords=[100, 200, 100, 200], time_range=[0,238],
#                        p0=[0.5, 45., 10., -1.1], plot=True, savefig=None)