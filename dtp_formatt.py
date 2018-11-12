"""
Farhad Allian
PhD Student
School of Mathematics & Statistics
Univesity of Sheffield (UK)

"""








import warnings  
warnings.filterwarnings("ignore")  

import matplotlib.pyplot as plt
import sunpy.map





"""

Finds the total number of FITS files in a given directory. 

Example call:
    
>> [in]: from dtp import userinput
>> [in]: userinput()
>> [out]: Where are your files stored?
>> [in]: Downloads
>> [out]: ~resulting filepath~
etc..


Note: Your files must be stored in the same directory as this script.

"""


import os

def listdir_fullpath(d):
        return [os.path.join(d, f) for f in os.listdir(d)]
    

location = 'D:/my_work/projects/Asymmetric_slab/Application_to_observation/Morton_2012_data/rosa_data'

file_path = os.path.abspath(location)

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
    
    print("\nTotal number of FITS files in this directory is %s. \n" %len(total_files))

    





def createmapcube(roi_x, roi_y, getmaps = False, normal_plot = False,  basic_plot=False,
                  draw_grid=False, vmin=0, vmax=2000, draw_limb=False,
                  submap=None, sobel_enhanced=False, get_limb=False):
    
    """ This function creates the data/map cube. 
    
    
    Inputs: 
        
        getmaps - whether you want to create the mapcube in the given directory.
        
        plot - creates the standard plot of the first image.
        
        basic_plot - same as above, but without axes and labels.
        
        draw_grid - plots the solar longitude and latitude and the limb.
        
        draw_limb - plots only the limb
        
        submap - whether you want to look at a specific region of interest. 
        
        vmin,vmax (optional) - intensity of images.
    
    Example call: 
        
        >> [in]: from dtp import createmapcube
        >> [in]: createmapcube(getmaps=True)
    """
    
    import sunpy
    from matplotlib import colors
    import matplotlib.pyplot as plt
    import astropy.units as u
    from skimage import exposure
    import numpy as npastroml
    import numpy as np
    from scipy import ndimage
#    from mgnn import mgn
    global header
    global total_maps, total_mapss,newmaps,newcubemaps
    header = []
    newmaps=[]
    newcubemaps=[]
    
    if getmaps == True:
        
        print("\nCreating datacube... ".format(len(total_files)))
        
        total_maps = sunpy.map.Map((str(file_path)+'/*ts'),sequence=True)
        
        
#        prepped=[]
#        for maps in total_maps:
#            
#            sunpy.instr.aia.aiaprep(maps)
        
#        total_maps = sunpy.map.Map(prepped,cube=True)
        
        
      #  total_maps = mapcube_solar_derotate(total_maps,clip=True)
        
#        for i in range(len(total_files)):
        
#            header.append(total_maps[i].meta)
            
#        for i in range(len(total_files)):
         
#           newmaps.append(mgn(total_maps[i],  gamma=3.2, g = 0.65 , h= 0.96) )
       
#        for i in range(len(total_files)):
#            
#            newcubemaps.append(sunpy.map.Map(newmaps[i,:,:], header[i]))
      
#        total_maps = sunpy.map.Map(newmaps,cube=True)
    
        
    if len(total_maps) == len(total_files):
        
        print("\nMapcube created.")

        
        if submap == True:
            
            print("\nCreating submap...")
            
            global sub_map
            sub_map= []
            
            from astropy.coordinates import SkyCoord
            import matplotlib.animation as animation
            
            
            global roi_x0, roi_y0, roi_x1, roi_y1
            roi_x0, roi_y0 = int(roi_x[0]), int(roi_y[0])
            roi_x1, roi_y1 = int(roi_x[1]), int(roi_y[1])
            
            global leftbox,rightbox
            leftbox = SkyCoord(roi_x0*u.arcsec, roi_y0*u.arcsec, frame = total_maps[0].coordinate_frame)

            rightbox = SkyCoord(roi_x1*u.arcsec, roi_y1*u.arcsec, frame = total_maps[0].coordinate_frame)
            
            print("\nCreating submap... ")
            
        
            newmaps=[]
            newcubemaps=[]
            
            #Mapcube has no attribute submap, so I'm doing it this way
            
#            for file in total_files:
#    
#                sub_map.append(sunpy.map.Map(file).submap(leftbox, rightbox))
#            
#            sub_map = sunpy.map.Map(sub_map, cube=True)
            
            for maps in total_maps:
                
                sub_map.append(maps.submap(leftbox,rightbox))
                
            sub_map = sunpy.map.Map(sub_map, cube=True)
        
            if len(total_files) == len(sub_map):
                
                print("\nPlotting region of interest...")
            
                sub_map[0].plot(norm=colors.Normalize(vmin=0,vmax=1))
                
                plt.title('ROI \n({0},{1}) ({2},{3}).'.format(roi_x0,roi_y0,roi_x1,roi_y1))
                
                plt.show()
                
                print("\nSubmap created.")
                
        if sobel_enhanced == True:
            
            global enhanced,enhancedd
            enhanced = []
            enhancedd =[]

            
            print("\nApplying Gaussian blur and Sobel enhancement to the submap... ")
            
            for file in total_files:
                
                m_smap = sunpy.map.Map(file).submap(leftbox, rightbox)
                
                header = m_smap.meta
                
                m_smap = exposure.equalize_hist(m_smap.data)
                
                m_smap = ndimage.gaussian_filter(m_smap,sigma=1.5)
            
                
                sx = ndimage.sobel(m_smap.data, axis=0, mode='constant')

                sy = ndimage.sobel(m_smap.data, axis=1, mode='constant')

                edge_enhanced_im = np.hypot(sx, sy)
                
              #  theta = np.arctan(sy/sx)

                enhanced.append(sunpy.map.Map(edge_enhanced_im, header))
                
          #      enhancedd.append(sunpy.map.Map(theta, m_smap.meta))
                
            
            enhanced= sunpy.map.Map(enhanced,cube=True)
       #     enhancedd = sunpy.map.Map(enhancedd,cube=True)
    
        
            print("\nEnhanced maps created.")
            
         #   enhanced[0].peek()
        
        if basic_plot == True:
            
            total_maps[0].peek(basic_plot=basic_plot, draw_limb=draw_limb,norm=colors.Normalize(vmin=vmin, vmax=vmax))
        
            print ("\nBasic plot created.")
        
        elif normal_plot == True:
            
            total_maps[0].peek(draw_grid=draw_grid, draw_limb=draw_limb,norm=colors.Normalize(vmin=vmin, vmax=vmax))
            
            print ("\nNormal plot created.")
    
    else:
        raise Exception("The number of files does not equal the number of maps. Please create mapcube again.")
            
createmapcube(roi_x=[0,100], roi_y=[0,100], getmaps = True, submap = False)














#def slices(submap=False):
#    
#    """ 
#    
#    This function prompts the user for their slice coordinates. 
#    
#    
#    It currently only works for one slice. Improve later.
#    
#    
#    Example call if you want a submap conversion from arcsec to pixels: 
#        
#        >> [in]: from dtp import slices
#        >> [in]: slices(submap=True)
#        >> [out]: "How many slices will you need? (integer)"
#        etc..
#  
#    """
#    
#    global reply,x0,x1,y0,y1,xp0,xp1,yp0,yp1
#    
#    reply = input("How many slices will you need? (integer) ")
#    
#    reply = int(reply)
#    
#
#    
#    #these functions convert arc seconds to pixels, for submap and total_maps
#    def a2p_sub(x, y):
#        
#            "Convert from arcseconds to pixels (submap)"
#        
#            xpixel = (x - sub_map[0].meta['crval1'])/sub_map[0].meta['cdelt1'] + (sub_map[0].meta['crpix1'] -1 )
#        
#            ypixel = (y - sub_map[0].meta['crval2'])/sub_map[0].meta['cdelt2'] + (sub_map[0].meta['crpix2'] -1 )
#
#            return xpixel,ypixel
#        
#    def a2p(x, y):
#        
#            "Convert from arcseconds to pixels (total_maps)"
#        
#            xpixel = (x - total_maps[0].meta['crval1'])/total_maps[0].meta['cdelt1'] + (total_maps[0].meta['crpix1'] -1 )
#        
#            ypixel = (y - total_maps[0].meta['crval2'])/total_maps[0].meta['cdelt2'] + (total_maps[0].meta['crpix2'] -1 )
#        
#            return xpixel,ypixel
#        
#        
#    
#    if reply == 1 and submap==True:
#        
#        print("\nCreating one slice...")
#        
#        
#        x0,y0 = (input("Enter your initial coordinates of your submap slice (arcsecs) in the form: x0,y0    ")).split(',')
#    
#        x1,y1 = (input('Enter your final coordinates of your submap slice (arcsec) in the form: x1,y1    ')).split(',')
#    
#        x0,y0,x1,y1 = int(x0), int(y0), int(x1), int(y1)
#        
#        print("""\nYour slice coordinates in arc seconds are:\n
#              (x0, y0) = ({0}, {1})''\n
#              (x1, y1) = ({2}, {3})''\n
#              """.format(x0,y0,x1,y1))
#        
#        print("\nConverting to pixels... ")
#            
#        
#        xp0,yp0 = a2p_sub(x0,y0)
#    
#        xp1,yp1 = a2p_sub(x1,y1)
#            
#        print("""\nDone. \n\nYour slice coordinates in pixels are:\n
#              (x0, y0) = ({0}, {1}) pix\n
#              (x1, y1) = ({2}, {3}) pix\n
#              """.format(xp0,yp0,xp1,yp1))
#
#    elif reply == 1 and submap == False:
#        
#        print("\nCreating one slice...")
#        
#        
#        x0,y0 = (input("Enter your initial coordinates of your slice (arcsecs) in the form: x0,y0    ")).split(',')
#    
#        x1,y1 = (input('Enter your final coordinates of your slice (arcsec) in the form: x1,y1    ')).split(',')
#    
#        x0,y0,x1,y1 = int(x0), int(y0), int(x1), int(y1)
#        
#        print("""\nYour slice coordinates in arc seconds are:\n
#              (x0, y0) = ({0}, {1})''\n
#              (x1, y1) = ({2}, {3})''\n
#              """.format(x0,y0,x1,y1))
#        
#        print("\nConverting to pixels... ")
#            
#        
#        xp0,yp0 = a2p(x0,y0)
#    
#        xp1,yp1 = a2p(x1,y1)
#            
#        print("""\nDone. \n\nYour slice coordinates in pixels are:\n
#              (x0, y0) = ({0}, {1}) pix\n
#              (x1, y1) = ({2}, {3}) pix\n
#              """.format(xp0,yp0,xp1,yp1))
#
#
#    #leave this for now
#    if reply > 1:
#        
#        global initial
#        initial= []
#        global final
#        final=[]
#        global inpixel
#        inpixel = []
#        global xvals
#        xvals = []
#        global yvals
#        yvals = []
#        global inpvals
#        inpvals=[]
#        global fipvals
#        fipvals=[]
#        global pixcoords
#        pixcoords=[]
#
#
#        print("\nCreating multiple slices...")
#        
#        ininput = input("Enter ALL of your initial coordinates in arcsecs in the form: x_i,y_i  (integer)  ").split(',')
#        initial.append(ininput)
#        fiinput = input("Enter ALL of your corresponding final coordinates in arcsecs in the form: x_f,y_f  (integer)  ").split(',')
#        final.append(fiinput)
#        
#        if 2*reply != len(initial[0]) and len(final[0]):
#            
#            raise Exception("\nPlease check you have entered the correct number of coordinates for your number of slices")
#            
#        print("\nConverting arc seconds to pixel coordinates...")
#        
#        for i in range(0, len(initial[0])):
#            initial[0][i] = int(initial[0][i])
#            final[0][i] = int(final[0][i])
#        
#           # inpvals = [x[0] for x in pixcoords]
#           # fipvals = [x[1] for x in pixcoords]
#        
#        #array slicing to get x and y values
#        xvals.append(initial[0][0:2*reply-1:reply])
#        xvals.append(final[0][0:2*reply-1:reply])
#        yvals.append(initial[0][1:2*reply:reply])
#        yvals.append(final[0][1:2*reply:reply])
#        
#    
#            
#        print("""\nYou wanted {0} slices.\n\nYour pixel coordinates are stored as a list of tuples.
#              """.format(reply))
#        
#slices(submap=False)










# slit coords in pixels
#slit_coords_x = [172, 172]
#slit_coords_y = [450, 499]

slit_coords_x = [611, 580]
slit_coords_y = [282, 371]

def distancetime(getintens = False, xfinal=slit_coords_x[0], xinitial=slit_coords_x[1], yfinal=slit_coords_y[0], yinitial=slit_coords_x[1]):
    
    """
    DISTAnces in pixels
    
    This function creates the distance-time plots along the user's slice coordiantes.
    
    getintens - whether you want to calculate the intensity values along the slice
    
    dt - whether you want to plot only distance time plot
    
    peaks, threhsold, box_size - finds the maximum intensity value on distance time plots. (For time series analysis)
    
    intensityslice - I = I(t,s(x,y))
    
    
    The zero point is at BOTTOM RIGHT  whereas the TOP LEFT is at 80 Mm.
    
    """
    import numpy as np
    import scipy.ndimage
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib import colors
    import warnings
#    from astroML.plotting import setup_text_plots
    from matplotlib.colors import LogNorm
    from scipy import ndimage
 
    
    
    warnings.filterwarnings("ignore") #this ignores some deprecated functions 
    
    if getintens == True:
    

        print("\nCalculating the intensity values along slit....")
        
        global intensity1
        intensity1 = []

        
        
        num = np.sqrt((xfinal - xinitial)**2 + (yfinal - yinitial)**2) #aia resolution relative. Calculated by doing space.max()/0.435
        
        global x,y
        x,y = np.linspace(xfinal,xinitial,num), np.linspace(yfinal,yinitial,num) #changed the order of xp1,xp0 and yp1,yp0

          
        for i in range(len(total_maps)):
    
            intensity1.append(scipy.ndimage.map_coordinates(np.transpose(total_maps[i].data),np.vstack((x,y))))
            
        cad = 7.68
        global time
        time = np.linspace(0, len(intensity1)*cad, len(intensity1))
              
        global space
        space = np.linspace(0, np.sqrt(((xfinal-xinitial)**2 + (yfinal-yinitial)**2))*50., len(intensity1[0])) # multiply by 0.428 instead to convert into Mm
        
        intensity1=np.array(intensity1)


        


    
    
    
    
 
    
import numpy as np
#time_window_pix = np.array([100,140]) #in pixels
time_window_pix = np.array([0,238]) #in pixels
time_window = time_window_pix * 7.68 # in s
time_frames = time_window_pix[1] - time_window_pix[0]
    
#       
distancetime(getintens=True)
#
plt.imshow(intensity1.T[:,time_window_pix[0]:time_window_pix[1]], aspect='auto', extent=[time_window[0], time_window[1], space.min(), space.max()], interpolation=None, origin='lower') #1100, 1500

plt.figure()
plt.imshow(total_maps[0].data, aspect='auto', interpolation=None, origin='lower')
plt.plot(slit_coords_x, slit_coords_y, color='white')




#################################
#Animate through time
import matplotlib.animation as animation

fig = plt.figure()

def updatefig(i):
    image_data = total_maps[i].data
    im.set_array(image_data)
    return im

im = plt.imshow(total_maps[time_window_pix[0]].data, aspect='auto', interpolation='gaussian', animated=True, origin='lower')

ani = animation.FuncAnimation(fig, updatefig, frames=time_frames, interval=200, repeat_delay = 1000) #interval = delay between each fram in ms
plt.plot(slit_coords_x, slit_coords_y, color='white')

