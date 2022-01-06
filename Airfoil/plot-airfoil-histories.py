#!/usr/bin/env python3
import numpy as np
import scipy.integrate as integrate
import scipy.fftpack   as fftpack
import os.path

import matplotlib
#matplotlib.use('pgf')
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams


# remove or set to True (default) to trigger exception
rc_xelatex = {'pgf.rcfonts': False} 
matplotlib.rcParams.update(rc_xelatex)

rc('text', usetex=True)
rc('font', family='serif')
rcParams.update({'figure.autolayout': True})


M1 = plt.figure(figsize=(12,4))
M1.set_facecolor('White')
ax_M1_1 = M1.add_subplot(131)
ax_M1_2 = M1.add_subplot(132)
ax_M1_3 = M1.add_subplot(133)

M2 = plt.figure(figsize=(12,4))
M2.set_facecolor('White')
ax_M2_1 = M2.add_subplot(131)
ax_M2_2 = M2.add_subplot(132)
ax_M2_3 = M2.add_subplot(133)



# AFRL Cases 
afrl_motions = ['M1', 'M2']

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
        #   e.g.
        #   nt     = 4
        #   nsteps = 3
        
            
        if ( isinstance(afrl_data, np.ndarray) ):
            # AFRL Data
            afrl_Fx   = afrl_data[:,0]
            afrl_Fy   = afrl_data[:,1]
            afrl_Fz   = afrl_data[:,2]
            afrl_Wint = afrl_data[:,3]
            afrl_mass = afrl_data[:,4]
            
            afrl_nx = len(afrl_Fy)
            xend    = 2.
            afrl_dx = 2./(afrl_nx-1)
            afrl_x  = np.linspace(0.,xend,afrl_nx)
            
            color = 'b--'
            if (motion == 'M1'):
                ax_M1_1.plot(afrl_x,afrl_Fx,  color,linewidth=1.0,label='AFRL')
                ax_M1_2.plot(afrl_x,afrl_Fy,  color,linewidth=1.0,label='AFRL')
                ax_M1_3.plot(afrl_x,afrl_Wint,color,linewidth=1.0,label='AFRL')
            elif (motion == 'M2'):
                ax_M2_1.plot(afrl_x,afrl_Fx,  color,linewidth=1.0,label='AFRL')
                ax_M2_2.plot(afrl_x,afrl_Fy,  color,linewidth=1.0,label='AFRL')
                ax_M2_3.plot(afrl_x,afrl_Wint,color,linewidth=1.0,label='AFRL')




# UM Cases 
UM_motions    = ['M1', 'M2']

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
        um_Fx   = um_data[:,4]
        um_Fy   = um_data[:,5]
        um_Wint = um_data[:,7]
        
        um_nx = len(um_Fy)
        xend    = 2.
        um_dx = 2./(um_nx-1)
        um_x  = np.linspace(0.,xend,um_nx)
        
    
        if (motion == 'M1'):
            ax_M1_1.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
            ax_M1_2.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
            ax_M1_3.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')

        elif (motion == 'M2'):
            ax_M2_1.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
            ax_M2_2.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
            ax_M2_3.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')



# ku cases 
ku_motions = ['M1', 'M2']

for imotion,motion in zip(range(len(ku_motions)),ku_motions):

    if motion == "M1":
        ku_file = 'KU/update_20211231/Re1000_Motion1.txt'

    elif motion == "M2":
        ku_file = 'KU/update_20211231/Re1000_Motion2.txt'


    # load data
    if (os.path.isfile(ku_file)):
        ku_data = np.loadtxt(ku_file,delimiter=',')
    else:
        print('ku data not found!')
        ku_data = false

    # reported data does not include initial time. 
    #   missing --- t1 --- t2 --- ... --- t_end
    #
    if ( isinstance(ku_data, np.ndarray) ):

        # insert info for t=0
        ku_data = np.append([[0., 0., 0., 0.]],ku_data,axis=0)

        # Remove duplicate entry for t_end or times greater than 2.
        if (ku_data[-1,0] == ku_data[-2,0]) or (ku_data[-1,0] > 2.00000001):
            print("KU: Removing duplicated final time")
            ku_data = np.delete(ku_data, -1, axis=0)

        # University of Kansas data
        ku_t    = ku_data[:,0]
        ku_fx   = ku_data[:,1]
        ku_fy   = ku_data[:,2]
        ku_wint = ku_data[:,3]
    
        if (motion == 'M1'):
            ax_M1_1.plot(ku_t,ku_fx,  'r--',linewidth=1.0,label='KU')
            ax_M1_2.plot(ku_t,ku_fy,  'r--',linewidth=1.0,label='KU')
            ax_M1_3.plot(ku_t,ku_wint,'r--',linewidth=1.0,label='KU')

        elif (motion == 'M2'):
            ax_M2_1.plot(ku_t,ku_fx,  'r--',linewidth=1.0,label='KU')
            ax_M2_2.plot(ku_t,ku_fy,  'r--',linewidth=1.0,label='KU')
            ax_M2_3.plot(ku_t,ku_wint,'r--',linewidth=1.0,label='KU')





# U. of Strasbourg cases 
us_motions = ['M1', 'M2']

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

        us_fx   = us_data[:,3]
        us_fy   = us_data[:,2]
        #us_wint = us_data[:,3]

        color = 'c--'
        if (motion == 'M1'):
            ax_M1_1.plot(us_t,us_fx,  color,linewidth=1.0,label='US')
            ax_M1_2.plot(us_t,us_fy,  color,linewidth=1.0,label='US')
            #ax_M1_3.plot(us_t,us_wint,color,linewidth=1.0,label='US')

        elif (motion == 'M2'):
            ax_M2_1.plot(us_t,us_fx,color,linewidth=1.0,label='US')
            ax_M2_2.plot(us_t,us_fy,color,linewidth=1.0,label='US')
            #ax_M2_3.plot(us_t,us_wint,color,linewidth=1.0,label='US')












#
#        if (group == 'ucb'):
#            if (case == 'case1'):
#                ucb_file = 'external/ucb/w1_2330020_Re10_case1_forces.dat'
#            elif (case == 'case2'):
#                ucb_file = 'no-file'
#            elif (case == 'case3'):
#                ucb_file = 'no-file'
#            elif (case == 'case4'):
#                ucb_file = 'no-file'
#
#
#            # Load data
#            if (os.path.isfile(ucb_file)):
#                ucb_data = np.loadtxt(ucb_file, skiprows=1)
#            else:
#                ucb_data = False
#
#            #print(type(ucb_data))
#            #print(ucb_data.shape)
#        
#            # Reported data does NOT include initial time, but does
#            # include final time
#            #   Missing --- t1 --- t2 --- t3
#            #
#            if ( isinstance(ucb_data, np.ndarray) ):
#
#                ucb_data = np.append([[0., 0., 0., 0., 0.]],ucb_data,axis=0)
#
#                # UC Berkeley Data
#                ucb_Fx   = ucb_data[:,1]
#                ucb_Fy   = ucb_data[:,2]
#                ucb_Wint = ucb_data[:,4]
#                
#                ucb_nx = len(ucb_Fy)
#                xend    = 2.
#                ucb_dx = 2./(ucb_nx-1)
#                ucb_x  = np.linspace(0.,xend,ucb_nx)
#
#                ucb_integrated_Fx = integrate.simps(ucb_Fx,   dx=ucb_dx)
#                ucb_integrated_Fy = integrate.simps(ucb_Fy,   dx=ucb_dx)
#                ucb_integrated_W  = integrate.simps(ucb_Wint, dx=ucb_dx)
#                #ucb_integrated_Fx = integrate.romb(ucb_Fx,   dx=ucb_dx)
#                #ucb_integrated_Fy = integrate.romb(ucb_Fy,   dx=ucb_dx)
#                #ucb_integrated_W  = integrate.romb(ucb_Wint, dx=ucb_dx)
#
#                if (motion == 'M1'):
#                    ax_c1_1.plot(ucb_x,ucb_Fx,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c1_2.plot(ucb_x,ucb_Fy,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c1_3.plot(ucb_x,ucb_Wint,'r-.',linewidth=1.0,label='UC-Berk')
#                elif (motion == 'M2'):
#                    ax_c2_1.plot(ucb_x,ucb_Fx,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c2_2.plot(ucb_x,ucb_Fy,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c2_3.plot(ucb_x,ucb_Wint,'r-.',linewidth=1.0,label='UC-Berk')
#                elif (motion == 'M3'):
#                    ax_c3_1.plot(ucb_x,ucb_Fx,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c3_2.plot(ucb_x,ucb_Fy,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c3_3.plot(ucb_x,ucb_Wint,'r-.',linewidth=1.0,label='UC-Berk')
#                elif (motion== 'M4'):
#                    ax_c4_1.plot(ucb_x,ucb_Fx,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c4_2.plot(ucb_x,ucb_Fy,  'r-.',linewidth=1.0,label='UC-Berk')
#                    ax_c4_3.plot(ucb_x,ucb_Wint,'r-.',linewidth=1.0,label='UC-Berk')
#
#    print("(Case, Re): ", " (", case, ",", re_number, ")")
#    print("Org, Fx, Fy, W, mass")
#    print("AFRL:", afrl_integrated_Fx, afrl_integrated_Fy, afrl_integrated_W, afrl_integrated_mass)
#    print("UM:",   um_integrated_Fx,   um_integrated_Fy,   um_integrated_W  )
#    print("UCB:",  ucb_integrated_Fx,  ucb_integrated_Fy,  ucb_integrated_W )


ax_M1_1.legend()
ax_M1_2.legend()
ax_M1_3.legend()

ax_M2_1.legend()
ax_M2_2.legend()
ax_M2_3.legend()


ax_M1_1.set_xlabel('Time')
ax_M1_2.set_xlabel('Time')
ax_M1_3.set_xlabel('Time')

ax_M2_1.set_xlabel('Time')
ax_M2_2.set_xlabel('Time')
ax_M2_3.set_xlabel('Time')

ax_M1_1.set_ylabel('Force-X (Euler)')
ax_M1_2.set_ylabel('Force-Y (Euler)')
ax_M1_3.set_ylabel('Power (Euler)')

ax_M2_1.set_ylabel('Force-X (Euler)')
ax_M2_2.set_ylabel('Force-Y (Euler)')
ax_M2_3.set_ylabel('Power (Euler)')


ax_M1_1.set_xlim((0.,2.))
ax_M1_2.set_xlim((0.,2.))
ax_M1_3.set_xlim((0.,2.))
#
ax_M2_1.set_xlim((0.,2.))
ax_M2_2.set_xlim((0.,2.))
ax_M2_3.set_xlim((0.,2.))

ax_M1_1.set_ylim((-0.5,0.5))
ax_M1_2.set_ylim((-2.5,2.5))
ax_M1_3.set_ylim((-2.5,2.5))
#
ax_M2_1.set_ylim((-2.5,4.5))
ax_M2_2.set_ylim((-4.5,4.5))
ax_M2_3.set_ylim((-6.5,1.5))


M1.savefig('Motion1_Airfoil_Histories.png', bbox_inches='tight', dpi=800)
M2.savefig('Motion2_Airfoil_Histories.png', bbox_inches='tight', dpi=800)

plt.show()




