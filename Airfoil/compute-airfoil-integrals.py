#!/usr/bin/env python3
import numpy as np
import scipy.integrate as integrate
import scipy.fftpack   as fftpack
import os.path

airfoil_motions = ['M1', 'M2']


# AFRL Cases 
afrl_motions    = ['M1', 'M2']
# Storage (motions,physics)
afrl_xI = np.zeros((2))
afrl_yI = np.zeros((2))
afrl_zI = np.zeros((2))
afrl_W  = np.zeros((2))
afrl_m  = np.zeros((2))


for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):

    motion1_truth = ['h3','p2','0.01']
    motion2_truth = ['h3','p2','0.01']

    if (motion == 'M1'):
        h_ref = motion1_truth[0]
        p_ref = motion1_truth[1]
        t_ref = motion1_truth[2]
    elif (motion == 'M2'):
        h_ref = motion2_truth[0]
        p_ref = motion2_truth[1]
        t_ref = motion2_truth[2]

    afrl_file = 'AFRL/'+motion+'-'+h_ref+'-'+p_ref+'/'+'airfoil-t'+t_ref

    # Load data
    if (os.path.isfile(afrl_file)):
        afrl_data = np.loadtxt(afrl_file, skiprows=1)
    else:
        afrl_data = False

    # Reported data includes initial and final time
    #   t0 --- t1 --- t2 --- t3
    #
    if ( isinstance(afrl_data, np.ndarray) ):
        # AFRL Data
        afrl_Fx   = afrl_data[:,0]
        afrl_Fy   = afrl_data[:,1]
        afrl_Fz   = afrl_data[:,2]
        afrl_Wint = afrl_data[:,3]
        afrl_mass = afrl_data[:,4]
        
        afrl_nx = len(afrl_Fy)
        xend    = 2.
        afrl_dx = 2./float(afrl_nx-1)
        afrl_x  = np.linspace(0.,xend,afrl_nx)
        
        afrl_xI[imotion] = integrate.simps(afrl_Fx,  dx=afrl_dx)
        afrl_yI[imotion] = integrate.simps(afrl_Fy,  dx=afrl_dx)
        afrl_zI[imotion] = integrate.simps(afrl_Fz,  dx=afrl_dx)
        afrl_W[ imotion] = integrate.simps(afrl_Wint,dx=afrl_dx)
        afrl_m[ imotion] = integrate.simps(afrl_mass-afrl_mass[0],dx=afrl_dx)









# UM Cases 
UM_motions    = ['M1', 'M2']
# Storage (motions,physics)
UM_xI = np.zeros((2))
UM_yI = np.zeros((2))
UM_zI = np.zeros((2))
UM_W  = np.zeros((2))
UM_m  = np.zeros((2))


for imotion,motion in zip(range(len(UM_motions)),UM_motions):

    um_file = 'UM/Tri-'+motion+'-''h3-p3/naca3_TimeHist.txt'

    # Load data
    if (os.path.isfile(um_file)):
        um_data = np.loadtxt(um_file, skiprows=37)
    else:
        print('UM data not found!')
        um_data = False

    # Reported data includes initial and final time
    #   t0 --- t1 --- t2 --- t3
    #
    #   nt     = 4
    #   nsteps = 3
    #   
    if ( isinstance(um_data, np.ndarray) ):
        # University of Michigan Data
        um_t    = um_data[:,1]
        um_Fx   = um_data[:,4]
        um_Fy   = um_data[:,5]
        um_Wint = um_data[:,7]

        
        UM_xI[imotion] = integrate.simps(um_Fx,  x=um_t)
        UM_yI[imotion] = integrate.simps(um_Fy,  x=um_t)
        UM_W[ imotion] = integrate.simps(um_Wint,x=um_t)



# KU Cases 
KU_motions = ['M1', 'M2']

# Storage (motions,physics)
KU_xI = np.zeros((2))
KU_yI = np.zeros((2))
KU_zI = np.zeros((2))
KU_W  = np.zeros((2))
KU_m  = np.zeros((2))


for imotion,motion in zip(range(len(KU_motions)),KU_motions):

    if motion == "M1":
        ku_file = 'KU/update_20211231/Re1000_Motion1.txt'

    elif motion == "M2":
        ku_file = 'KU/update_20211231/Re1000_Motion2.txt'


    # Load data
    if (os.path.isfile(ku_file)):
        ku_data = np.loadtxt(ku_file,delimiter=',')
    else:
        print('KU data not found!')
        ku_data = False

    # Reported data does NOT include initial time. 
    #   Missing --- t1 --- t2 --- ... --- t_end
    #
    if ( isinstance(ku_data, np.ndarray) ):

        # Insert info for t=0
        ku_data = np.append([[0., 0., 0., 0.]],ku_data,axis=0)

        # Remove duplicate entry for t_end or times greater than 2.
        if (ku_data[-1,0] == ku_data[-2,0]) or (ku_data[-1,0] > 2.00000001):
            print("KU: Removing duplicated final time")
            ku_data = np.delete(ku_data, -1, axis=0)

        # University of Kansas Data
        ku_t    = ku_data[:,0]
        ku_Fx   = ku_data[:,1]
        ku_Fy   = ku_data[:,2]
        ku_Wint = ku_data[:,3]
        
        KU_xI[imotion] = integrate.simps(ku_Fx,  x=ku_t)
        KU_yI[imotion] = integrate.simps(ku_Fy,  x=ku_t)
        KU_W[ imotion] = integrate.simps(ku_Wint,x=ku_t)



# U. of Strasbourg cases 
us_motions = ['M1', 'M2']

# Storage (motions,physics)
US_xI = np.zeros((2))
US_yI = np.zeros((2))
US_zI = np.zeros((2))
US_W  = np.zeros((2))
US_m  = np.zeros((2))

for imotion,motion in zip(range(len(us_motions)),us_motions):

    if motion == "M1":
        us_file = 'US/forces_T_DB3.dat'

    elif motion == "M2":
        us_file = 'US/forces_TR_DB3.dat'


    # load data
    if (os.path.isfile(us_file)):
        us_data = np.loadtxt(us_file)
    else:
        print('us data not found!')
        us_data = false

    # reported data does not include initial time. 
    #   missing --- t1 --- t2 --- ... --- t_end
    #
    if ( isinstance(us_data, np.ndarray) ):

        # insert info for t=0
        us_data = np.append([[0., 0., 0., 0., 0.]],us_data,axis=0)

        # university of kansas data
        us_t    = us_data[:,1]

        us_Fx   = us_data[:,3]
        us_Fy   = us_data[:,2]
        #us_wint = us_data[:,3]


        US_xI[imotion] = integrate.simps(us_Fx,  x=us_t)
        US_yI[imotion] = integrate.simps(us_Fy,  x=us_t)
        #US_W[ imotion] = integrate.simps(us_Wint,x=us_t)







# UCB Cases 
UCB_motions    = ['M1', 'M2']
# Storage (motions,physics)
UCB_xI = np.zeros((2))
UCB_yI = np.zeros((2))
UCB_zI = np.zeros((2))
UCB_W  = np.zeros((2))
UCB_m  = np.zeros((2))


for imotion,motion in zip(range(len(UCB_motions)),UCB_motions):

    ucb_file = 'UCB/'+motion+'_ref2p3.dat'

    # Load data
    if (os.path.isfile(ucb_file)):
        ucb_data = np.loadtxt(ucb_file, skiprows=1)
    else:
        print(ucb_file)
        print('UCB data not found!')
        ucb_data = False



    # Reported data includes initial and final time
    #   t0 --- t1 --- t2 --- t3
    #
    #   nt     = 4
    #   nsteps = 3
    #   
    if ( isinstance(ucb_data, np.ndarray) ):
        # U.C. Berkeley Data
        ucb_t    = ucb_data[:,0]
        ucb_Fx   = ucb_data[:,1]
        ucb_Fy   = ucb_data[:,2]
        ucb_Wint = ucb_data[:,3]

        UCB_xI[imotion] = integrate.simps(ucb_Fx,  x=ucb_t)
        UCB_yI[imotion] = integrate.simps(ucb_Fy,  x=ucb_t)
        UCB_W[ imotion] = integrate.simps(ucb_Wint,x=ucb_t)







print(" ")
print("Org.  ", "X-Impulse", "Y-Impulse", "Work")
for imotion,motion in zip(range(len(airfoil_motions)),airfoil_motions):
    print("---------------------------------------------------------------- ")
    print(" ")
    print(motion)
    print("AFRL:  ", afrl_xI[imotion], afrl_yI[imotion], afrl_W[imotion])
    print("UMich: ",   UM_xI[imotion],   UM_yI[imotion],   UM_W[imotion])
    print("KU:    ",   KU_xI[imotion],   KU_yI[imotion],   KU_W[imotion])
    print("UCB:   ",  UCB_xI[imotion],  UCB_yI[imotion],  UCB_W[imotion]) 
    print("US:    ",   US_xI[imotion],   US_yI[imotion],  "Not reported") 




