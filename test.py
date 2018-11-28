# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield
"""

import time_distance as td

file_path = 'D:/my_work/projects/Asymmetric_slab/Application_to_observation/Morton_2012_data/rosa_data'

morton12 = td.Full_map(file_path)

#morton12.animate(slit_coords=[100, 200, 100, 200])
#morton12.total_maps[0].peek()
d100 = morton12.total_maps[100].data
d0 = d100 = morton12.total_maps[0].data