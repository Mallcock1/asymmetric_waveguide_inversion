# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield

Fibril boundary data from Morton et al 2012 ROSA data
"""

import numpy as np

fibril_1 = False
fibril_2 = False
fibril_3 = False
fibril_4 = False

fibril_1 = True
#fibril_2 = True
#fibril_3 = True
#fibril_4 = True

############################################################
# Fibril 1 (pixel counting from Fig 4 in Morton et al 2012)

# function to stretch a given slit by a factor
def stretch_slit(slit_coords, factor):
    x0 = slit_coords[0]
    x1 = slit_coords[1]
    y0 = slit_coords[2]
    y1 = slit_coords[3]

    xmid = (x1 + x0) / 2
    ymid = (y1 + y0) / 2

    x0new = (x0 - xmid)*factor + xmid
    x1new = (x1 - xmid)*factor + xmid
    y0new = (y0 - ymid)*factor + ymid
    y1new = (y1 - ymid)*factor + ymid
    return [x0new, x1new, y0new, y1new]


if fibril_1 is True:
    print("You have chosen fibril_1")
    fibril_number = "fibril1"

    # width_pix and yb_pix from pixel counting.
    # Width of fibril in pixels
    width_pix = np.array([16.25,14,11,7.5,8,7,7,5.75,7.5,8.5,9.25,20.25,11.5,11.5,
                          11.75,9.5,8.5,7.5,7.25,7.25,6.25,5.5,7,6.5,7.75,8,6.75,6,
                          7.75,9,12,12.5,11.5,9])

    # y-distance from bottom of domain in pixels
    yb_pix = np.array([24.75,27.25,29.5,31.5,31,31.75,32.25,33.75,34,34.5,34.75,30,
                       35.5,36,36.25,37,36.75,36.5,37,37,37.75,38.25,38,39.5,41,
                       42.5,43.5,44,43,41.25,39,38,37,36.75])
    
#    # width_pix and yb_pix from data.
#    # Width of fibril in pixels
#    width_pix = np.array([])
#
#    # y-distance from bottom of domain in pixels
#    yb_pix = np.array([])

    yt_pix = yb_pix + width_pix

    unit = "pix"

    smooth_indices = [11]

    slit_coords = [420, 424, 615.6, 565.6]  # This is what the boundary data above is from: [425, 429, 616, 566]
    
#    slit_coords = stretch_slit(slit_coords, 1.5)
    
    time_range = slice(130, 225)  # slice(151, 225)  # in frames

    trend_range = slice(0, 22)

    cad = 7.68  # temporal resolution (cadence) of the ROSA instrument in s
    t_start = 1000 + 19*cad

    t_vals = np.arange(len(yb_pix))*cad*2 + t_start

    # Degree of the trend polynomial
    N = 2
    p0 = [100., 0.01, 0.]

    vA_guess = 30.  # 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.2  # R1 := rho_1 / rho_0
    R2 = 0.1  # R2 := rho_2 / rho_0
    c_phase = 47  # Morton found this: 71.
    mode = "saus"


###########################################################
# Fibril 2

if fibril_2 is True:
    print("You have chosen fibril_2")
    fibril_number = "fibril2"

    yb_pix = np.array([31.924199989670257,
                       33.55552590115389,
                       34.41517726955584,
                       34.7120289190419,
                       34.23829906223432,
                       36.66662332698089,
                       35.2013053417571,
                       34.18845386386369,
                       32.15205006253402,
                       30.948884235841664,
                       30.015207664499478,
                       31.666625618477774,
                       33.26768346386292,
                       34.667103455342,
                       35.65406299508912,
                       36.627021326950846,
                       36.07314345317778,
                       35.57445067302206,
                       34.89150223159021,
                       35.82574242720006,
                       36.64066157798094,
                       35.23837102073869,
                       34.37018316424221,
                       32.664718718103394,
                       21.55025711496983,
                       35.04824155568525,
                       35.94151809394504,
                       36.71602354728906,
                       35.68302817032468,
                       32.678981446207715,
                       29.825078594481887,
                       25.806759187549858,
                       25.83975342953978,
                       25.42683607728384,
                       28.385598098903596,
                       36.42959480492649,
                       40.70615699649271,
                       42.518128026445254,
                       42.78616759209693,
                       42.820821053781074,
                       41.28432891282563,
                       42.829484860958395,
                       43.398429713266225,
                       44.19912878182627,
                       44.0397485130752,
                       44.87388405806714,
                       45.22014434694579,
                       42.82658875782191,
                       41.884265838756384,
                       42.446205668842715,
                       35.08171411216903,
                       34.64563762830989,
                       31.974247530139184,
                       30.262730911052603,
                       23.618179512382202,
                       24.18149285504029,
                       23.05146106608389,
                       23.70969993466525,
                       25.116984205488933,
                       26.244308025355917,
                       27.62697563186082,
                       29.034941060109652,
                       30.94132433944659,
                       31.608658385284294,
                       34.51510732768399,
                       34.84517765809254,
                       37.67278324879001,
                       40.734339206613306,
                       43.93075142544651,
                       48.90253611251976,
                       49.00280796323902,
                       46.211029628742686])

    yt_pix = np.array([51.80763726015182,
                       50.717543662763696,
                       51.30231378375507,
                       51.13195059182715,
                       50.7695989016585,
                       50.98978454292617,
                       50.054207188130924,
                       50.89078865891877,
                       49.70558113451063,
                       51.24337632954297,
                       50.40812745519798,
                       51.077559659228264,
                       52.496292501052366,
                       53.37516427568582,
                       52.507642617759885,
                       52.94848439419122,
                       52.22905622499805,
                       51.53296975578987,
                       50.86212627539905,
                       49.58894950452258,
                       49.95455175748497,
                       49.05788859289485,
                       49.25524482708801,
                       48.981529111560654,
                       51.97447266699731,
                       49.851830948916785,
                       49.80669156909702,
                       50.90121338571256,
                       50.72658825524942,
                       52.873124122416314,
                       55.36860368285725,
                       57.72146466967121,
                       55.69513398421209,
                       56.14157246397944,
                       56.56495147830849,
                       55.49677892309504,
                       52.78278832783167,
                       52.854646294536764,
                       53.42137193432748,
                       53.624193625896964,
                       52.847002046642004,
                       52.87870729985459,
                       53.58093596756224,
                       54.19164866203635,
                       53.93340937430367,
                       54.7161591166006,
                       55.78697972830223,
                       55.780068628758585,
                       56.31523016018641,
                       58.48096677469393,
                       61.349133798675915,
                       62.32131624880845,
                       62.82660738719072,
                       63.878776465000996,
                       66.4589966314867,
                       65.20785476450658,
                       64.34961709327622,
                       60.42092824853303,
                       59.13869901354225,
                       60.083655691446864,
                       61.970265512835184,
                       61.34349376876472,
                       60.58423296749396,
                       60.46826789910539,
                       60.66497745363819,
                       59.19529214392901,
                       59.58094151195334,
                       58.28726509215321,
                       59.474063855227264,
                       59.95751548084363,
                       60.539867522023556,
                       60.711113499861504])

    unit = "pix"

    smooth_indices = [24]

    slit_coords = [591, 580, 265, 371]

    time_range = slice(0, 71)
    trend_range = slice(25, -1)

    cad = 7.68  # temporal resolution (cadence) of the ROSA instrument in s
    t_start = 0

    t_vals = np.arange(len(yb_pix))*cad + t_start

    # Degree of the trend polynomial
    N = 2
    p0 = [100., 0.01, 0.]

    vA_guess = 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.2  # R1 := rho_1 / rho_0
    R2 = 0.1  # R2 := rho_2 / rho_0
    c_phase = 71.
    mode = "saus"


###########################################################
# Fibril 3

if fibril_3 is True:
    print("You have chosen fibril_3")
    fibril_number = "fibril3"

    yb = np.array([743.2236682664131,
                   991.8614522308072,
                   881.8147596969167,
                   1056.9047351777654,
                   1134.7313046032928,
                   987.7467736626846,
                   1151.5202730595072,
                   1189.835318226618,
                   1045.1070381110846,
                   1147.7058796758395,
                   1377.5524712238987,
                   1440.4657413538885,
                   1418.5870285354708,
                   1253.10542296539,
                   1235.3747408615347,
                   1157.7911499880506,
                   1094.8720022910343,
                   1084.953726669772,
                   1109.5478614192918,
                   1073.8630650931448,
                   1055.0850857946284,
                   1039.2006300869878,
                   1071.481554068373,
                   981.9051279796081,
                   829.366471095187,
                   790.0970757685475,
                   871.3840204065889,
                   822.4022679562662,
                   815.7178906295618,
                   816.9326312208963,
                   769.7179272822918,
                   820.0331120599675,
                   843.8076265796287,
                   835.090627160881,
                   830.3280508600757,
                   784.5851778322576,
                   822.2398619507769,
                   828.9004214118993,
                   843.5458384839708,
                   804.7532231897611,
                   734.0761235593882,
                   887.0619108983803,
                   883.4911899445658])

    yt = np.array([2837.133940867595,
                   2775.5618890964174,
                   2758.0127280157667,
                   2737.7897138702756,
                   2789.484694663397,
                   2596.6094064268345,
                   2336.5390585753685,
                   2337.4323530221664,
                   2162.4998706300785,
                   2324.1565480271884,
                   2578.542590526814,
                   2690.9867649356383,
                   2712.550790427497,
                   2707.1927777796905,
                   2627.638288623935,
                   2687.5531837896624,
                   2448.0370551244255,
                   2277.1181927632606,
                   2121.6067633432685,
                   2145.409364637937,
                   1957.320558714092,
                   1963.6931435291406,
                   1750.644079420869,
                   1815.3579322938645,
                   1664.6849313230562,
                   1833.8474002695275,
                   1869.1652881884006,
                   1835.559215824414,
                   1921.9307475742746,
                   1720.6701563132392,
                   1546.8136377005528,
                   1450.466786469136,
                   1500.7918520500343,
                   1526.6388356829407,
                   1616.3791593470662,
                   1872.4126146509616,
                   1860.5447030306777,
                   1776.596882778248,
                   1775.1215258447176,
                   1699.6069573593586,
                   1700.084686363544,
                   1594.3533726967835,
                   1617.6993877802988])

    unit = "km"

    smooth_indices = [24]

    slit_coords = [775, 775, 430, 490]

    #time_window_frames = np.array([100,140]) #in frames 
    time_window_frames = np.array([132,175]) #in frames

    time_range = slice(10, -1)

    cad = 7.68  # temporal resolution (cadence) of the ROSA instrument in s
    t_start = 0

    t_vals = np.arange(len(yb_pix))*cad + t_start

    # Degree of the trend polynomial
    N = 2
    p0 = [100., 0.01, 0.]

    vA_guess = 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.2  # R1 := rho_1 / rho_0
    R2 = 0.1  # R2 := rho_2 / rho_0
    c_phase = 71.
    mode = "saus"
