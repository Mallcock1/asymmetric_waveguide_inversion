# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 14:05:52 2018

@author: Matt
"""

>>> Writer = animation.writers['ffmpeg']   
>>> writer = Writer(fps=10, metadata=dict(artist='SunPy'), bitrate=1800)   
>>>
>>> ani.save('mapcube_animation.mp4', writer=writer) 