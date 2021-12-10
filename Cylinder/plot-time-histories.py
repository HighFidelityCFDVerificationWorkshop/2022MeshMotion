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


M1 = plt.figure(figsize=(10,9))
M1.set_facecolor('White')
ax_M1_1 = M1.add_subplot(331)
ax_M1_2 = M1.add_subplot(332)
ax_M1_3 = M1.add_subplot(333)
ax_M1_4 = M1.add_subplot(334)
ax_M1_5 = M1.add_subplot(335)
ax_M1_6 = M1.add_subplot(336)
ax_M1_7 = M1.add_subplot(337)
ax_M1_8 = M1.add_subplot(338)
ax_M1_9 = M1.add_subplot(339)

M2 = plt.figure(figsize=(10,9))
M2.set_facecolor('White')
ax_M2_1 = M2.add_subplot(331)
ax_M2_2 = M2.add_subplot(332)
ax_M2_3 = M2.add_subplot(333)
ax_M2_4 = M2.add_subplot(334)
ax_M2_5 = M2.add_subplot(335)
ax_M2_6 = M2.add_subplot(336)
ax_M2_7 = M2.add_subplot(337)
ax_M2_8 = M2.add_subplot(338)
ax_M2_9 = M2.add_subplot(339)

M3 = plt.figure(figsize=(10,9))
M3.set_facecolor('White')
ax_M3_1 = M3.add_subplot(331)
ax_M3_2 = M3.add_subplot(332)
ax_M3_3 = M3.add_subplot(333)
ax_M3_4 = M3.add_subplot(334)
ax_M3_5 = M3.add_subplot(335)
ax_M3_6 = M3.add_subplot(336)
ax_M3_7 = M3.add_subplot(337)
ax_M3_8 = M3.add_subplot(338)
ax_M3_9 = M3.add_subplot(339)

M4 = plt.figure(figsize=(10,9))
M4.set_facecolor('White')
ax_M4_1 = M4.add_subplot(331)
ax_M4_2 = M4.add_subplot(332)
ax_M4_3 = M4.add_subplot(333)
ax_M4_4 = M4.add_subplot(334)
ax_M4_5 = M4.add_subplot(335)
ax_M4_6 = M4.add_subplot(336)
ax_M4_7 = M4.add_subplot(337)
ax_M4_8 = M4.add_subplot(338)
ax_M4_9 = M4.add_subplot(339)



# AFRL Cases 
afrl_grids   = ['h0', 'h1', 'h2']
afrl_motions = ['M1', 'M2', 'M3', 'M4']
afrl_physics = ['Re10','Re1000']
afrl_orders  = ['p1', 'p2', 'p3']
afrl_times   = ['0.1', '0.01', '0.001']

for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):


        if physic == 'Re10':
            motion1_truth = ['h0','p3','0.001']
            motion2_truth = ['h0','p3','0.001']
            motion3_truth = ['h0','p3','0.001']
            motion4_truth = ['h0','p3','0.001']
        elif physic == 'Re1000':
            motion1_truth = ['h0','p3','0.001']
            motion2_truth = ['h0','p3','0.001']
            motion3_truth = ['h0','p3','0.01']
            motion4_truth = ['h0','p2','0.01']
        else:
            print("ERROR!!!!!")


        if (motion == 'M1'):
            h_ref = motion1_truth[0]
            p_ref = motion1_truth[1]
            t_ref = motion1_truth[2]
        elif (motion == 'M2'):
            h_ref = motion2_truth[0]
            p_ref = motion2_truth[1]
            t_ref = motion2_truth[2]
        elif (motion == 'M3'):
            h_ref = motion3_truth[0]
            p_ref = motion3_truth[1]
            t_ref = motion3_truth[2]
        elif (motion == 'M4'):
            h_ref = motion4_truth[0]
            p_ref = motion4_truth[1]
            t_ref = motion4_truth[2]

        afrl_file = 'AFRL/'+motion+'-'+physic+'-'+h_ref+'-'+p_ref+'/'+'cyl-t'+t_ref
    
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
            
            if (motion == 'M1'):
                if (physic == 'ReInf'):
                    ax_M1_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M1_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M1_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re1000'):
                    ax_M1_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M1_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M1_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re10'):
                    ax_M1_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M1_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M1_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
            elif (motion == 'M2'):
                if (physic == 'ReInf'):
                    ax_M2_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M2_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M2_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re1000'):
                    ax_M2_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M2_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M2_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re10'):
                    ax_M2_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M2_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M2_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
            elif (motion == 'M3'):
                if (physic == 'ReInf'):
                    ax_M3_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M3_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M3_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re1000'):
                    ax_M3_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M3_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M3_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re10'):
                    ax_M3_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M3_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M3_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
            elif (motion == 'M4'):
                if (physic == 'ReInf'):
                    ax_M4_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M4_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M4_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re1000'):
                    ax_M4_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M4_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M4_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
                elif (physic == 'Re10'):
                    ax_M4_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
                    ax_M4_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
                    ax_M4_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')





# UM Cases 
UM_motions = ['M1', 'M2', 'M3', 'M4']
UM_physics = ['Re1','Re2','Re3']
UM_physics_re = ['Re1000','Re10','ReInf']

for imotion,motion in zip(range(len(UM_motions)),UM_motions):
    for iphysic,physic in zip(range(len(UM_physics)),UM_physics):

        um_file = 'UM/Quad-'+motion+'-'+physic+'-h2-p4/cyl3_TimeHist.txt'

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
                if (physic == 'Re3'):
                    ax_M1_1.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M1_2.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M1_3.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re2'):
                    ax_M1_4.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M1_5.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M1_6.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re1'):
                    ax_M1_7.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M1_8.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M1_9.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')

            elif (motion == 'M2'):
                if (physic == 'Re3'):
                    ax_M2_1.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M2_2.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M2_3.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re2'):
                    ax_M2_4.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M2_5.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M2_6.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re1'):
                    ax_M2_7.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M2_8.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M2_9.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')

            elif (motion == 'M3'):
                if (physic == 'Re3'):
                    ax_M3_1.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M3_2.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M3_3.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re2'):
                    ax_M3_4.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M3_5.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M3_6.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re1'):
                    ax_M3_7.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M3_8.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M3_9.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')

            elif (motion == 'M4'):
                if (physic == 'Re3'):
                    ax_M4_1.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M4_2.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M4_3.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re2'):
                    ax_M4_4.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M4_5.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M4_6.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')
                elif (physic == 'Re1'):
                    ax_M4_7.plot(um_x,um_Fx,  'y--',linewidth=1.0,label='UMich')
                    ax_M4_8.plot(um_x,um_Fy,  'y--',linewidth=1.0,label='UMich')
                    ax_M4_9.plot(um_x,um_Wint,'y--',linewidth=1.0,label='UMich')


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
#                # University of Michigan Data
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
ax_M1_4.legend()
ax_M1_5.legend()
ax_M1_6.legend()
ax_M1_7.legend()
ax_M1_8.legend()
ax_M1_9.legend()

ax_M2_1.legend()
ax_M2_2.legend()
ax_M2_3.legend()
ax_M2_4.legend()
ax_M2_5.legend()
ax_M2_6.legend()
ax_M2_7.legend()
ax_M2_8.legend()
ax_M2_9.legend()

ax_M3_1.legend()
ax_M3_2.legend()
ax_M3_3.legend()
ax_M3_4.legend()
ax_M3_5.legend()
ax_M3_6.legend()
ax_M3_7.legend()
ax_M3_8.legend()
ax_M3_9.legend()

ax_M4_1.legend()
ax_M4_2.legend()
ax_M4_3.legend()
ax_M4_4.legend()
ax_M4_5.legend()
ax_M4_6.legend()
ax_M4_7.legend()
ax_M4_8.legend()
ax_M4_9.legend()

ax_M1_1.set_xlabel('Time')
ax_M1_2.set_xlabel('Time')
ax_M1_3.set_xlabel('Time')
ax_M1_4.set_xlabel('Time')
ax_M1_5.set_xlabel('Time')
ax_M1_6.set_xlabel('Time')
ax_M1_7.set_xlabel('Time')
ax_M1_8.set_xlabel('Time')
ax_M1_9.set_xlabel('Time')

ax_M2_1.set_xlabel('Time')
ax_M2_2.set_xlabel('Time')
ax_M2_3.set_xlabel('Time')
ax_M2_4.set_xlabel('Time')
ax_M2_5.set_xlabel('Time')
ax_M2_6.set_xlabel('Time')
ax_M2_7.set_xlabel('Time')
ax_M2_8.set_xlabel('Time')
ax_M2_9.set_xlabel('Time')

ax_M3_1.set_xlabel('Time')
ax_M3_2.set_xlabel('Time')
ax_M3_3.set_xlabel('Time')
ax_M3_4.set_xlabel('Time')
ax_M3_5.set_xlabel('Time')
ax_M3_6.set_xlabel('Time')
ax_M3_7.set_xlabel('Time')
ax_M3_8.set_xlabel('Time')
ax_M3_9.set_xlabel('Time')

ax_M4_1.set_xlabel('Time')
ax_M4_2.set_xlabel('Time')
ax_M4_3.set_xlabel('Time')
ax_M4_4.set_xlabel('Time')
ax_M4_5.set_xlabel('Time')
ax_M4_6.set_xlabel('Time')
ax_M4_7.set_xlabel('Time')
ax_M4_8.set_xlabel('Time')
ax_M4_9.set_xlabel('Time')

ax_M1_1.set_ylabel('Force-X')
ax_M1_2.set_ylabel('Force-Y')
ax_M1_3.set_ylabel('Power')
ax_M1_4.set_ylabel('Force-X')
ax_M1_5.set_ylabel('Force-Y')
ax_M1_6.set_ylabel('Power')
ax_M1_7.set_ylabel('Force-X')
ax_M1_8.set_ylabel('Force-Y')
ax_M1_9.set_ylabel('Power')

ax_M2_1.set_ylabel('Force-X')
ax_M2_2.set_ylabel('Force-Y')
ax_M2_3.set_ylabel('Power')
ax_M2_4.set_ylabel('Force-X')
ax_M2_5.set_ylabel('Force-Y')
ax_M2_6.set_ylabel('Power')
ax_M2_7.set_ylabel('Force-X')
ax_M2_8.set_ylabel('Force-Y')
ax_M2_9.set_ylabel('Power')

ax_M3_1.set_ylabel('Force-X')
ax_M3_2.set_ylabel('Force-Y')
ax_M3_3.set_ylabel('Power')
ax_M3_4.set_ylabel('Force-X')
ax_M3_5.set_ylabel('Force-Y')
ax_M3_6.set_ylabel('Power')
ax_M3_7.set_ylabel('Force-X')
ax_M3_8.set_ylabel('Force-Y')
ax_M3_9.set_ylabel('Power')

ax_M4_1.set_ylabel('Force-X')
ax_M4_2.set_ylabel('Force-Y')
ax_M4_3.set_ylabel('Power')
ax_M4_4.set_ylabel('Force-X')
ax_M4_5.set_ylabel('Force-Y')
ax_M4_6.set_ylabel('Power')
ax_M4_7.set_ylabel('Force-X')
ax_M4_8.set_ylabel('Force-Y')
ax_M4_9.set_ylabel('Power')



ax_M1_1.set_xlim((0.,2.))
ax_M1_2.set_xlim((0.,2.))
ax_M1_3.set_xlim((0.,2.))
ax_M1_4.set_xlim((0.,2.))
ax_M1_5.set_xlim((0.,2.))
ax_M1_6.set_xlim((0.,2.))
ax_M1_7.set_xlim((0.,2.))
ax_M1_8.set_xlim((0.,2.))
ax_M1_9.set_xlim((0.,2.))
#
ax_M2_1.set_xlim((0.,2.))
ax_M2_2.set_xlim((0.,2.))
ax_M2_3.set_xlim((0.,2.))
ax_M2_4.set_xlim((0.,2.))
ax_M2_5.set_xlim((0.,2.))
ax_M2_6.set_xlim((0.,2.))
ax_M2_7.set_xlim((0.,2.))
ax_M2_8.set_xlim((0.,2.))
ax_M2_9.set_xlim((0.,2.))
#
ax_M3_1.set_xlim((0.,2.))
ax_M3_2.set_xlim((0.,2.))
ax_M3_3.set_xlim((0.,2.))
ax_M3_4.set_xlim((0.,2.))
ax_M3_5.set_xlim((0.,2.))
ax_M3_6.set_xlim((0.,2.))
ax_M3_7.set_xlim((0.,2.))
ax_M3_8.set_xlim((0.,2.))
ax_M3_9.set_xlim((0.,2.))
#
ax_M4_1.set_xlim((0.,2.))
ax_M4_2.set_xlim((0.,2.))
ax_M4_3.set_xlim((0.,2.))
ax_M4_4.set_xlim((0.,2.))
ax_M4_5.set_xlim((0.,2.))
ax_M4_6.set_xlim((0.,2.))
ax_M4_7.set_xlim((0.,2.))
ax_M4_8.set_xlim((0.,2.))
ax_M4_9.set_xlim((0.,2.))

#ax_M1_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax_M1_1.set_ylim((-0.1e-03,0.1e-03))
ax_M1_2.set_ylim((-1.0,2.5))
ax_M1_3.set_ylim((-1.0,1.0))
ax_M1_4.set_ylim((-0.1e-03,0.1e-03))
ax_M1_5.set_ylim((-1.0,2.5))
ax_M1_6.set_ylim((-1.0,1.0))
ax_M1_7.set_ylim((-0.1e-03,0.1e-03))
ax_M1_8.set_ylim((-1.0,2.5))
ax_M1_9.set_ylim((-1.0,1.0))
#
ax_M2_1.set_ylim((-0.1e-03,0.1e-03))
ax_M2_2.set_ylim((-0.1e-03,0.1e-03))
ax_M2_3.set_ylim((-1.0,1.0))
ax_M2_4.set_ylim((-0.1e-03,0.1e-03))
ax_M2_5.set_ylim((-0.1e-03,0.1e-03))
ax_M2_6.set_ylim((-1.0,1.0))
ax_M2_7.set_ylim((-0.1e-03,0.1e-03))
ax_M2_8.set_ylim((-0.1e-03,0.1e-03))
ax_M2_9.set_ylim((-1.0,1.0))
#
ax_M3_1.set_ylim((-0.1e-03,0.1e-03))
ax_M3_2.set_ylim((-0.1e-03,0.1e-03))
ax_M3_3.set_ylim((-0.05,0.05))
ax_M3_4.set_ylim((-0.1e-03,0.1e-03))
ax_M3_5.set_ylim((-0.1e-03,0.1e-03))
ax_M3_6.set_ylim((-0.05,0.05))
ax_M3_7.set_ylim((-0.1e-03,0.1e-03))
ax_M3_8.set_ylim((-0.1e-03,0.1e-03))
ax_M3_9.set_ylim((-0.05,0.05))
#
ax_M4_1.set_ylim((-0.1,0.1))
ax_M4_2.set_ylim((-1.0,2.5))
ax_M4_3.set_ylim((-1.5,2.0))
ax_M4_4.set_ylim((-0.1,0.1))
ax_M4_5.set_ylim((-1.0,2.5))
ax_M4_6.set_ylim((-1.5,2.0))
ax_M4_7.set_ylim((-0.1,0.1))
ax_M4_8.set_ylim((-1.0,2.5))
ax_M4_9.set_ylim((-1.5,2.0))



#case1.savefig('case1_re10.png', bbox_inches='tight', dpi=800)
#case2.savefig('case2_re10.png', bbox_inches='tight', dpi=800)
#case3.savefig('case3_re10.png', bbox_inches='tight', dpi=800)
#case4.savefig('case4_re10.png', bbox_inches='tight', dpi=800)

plt.show()




