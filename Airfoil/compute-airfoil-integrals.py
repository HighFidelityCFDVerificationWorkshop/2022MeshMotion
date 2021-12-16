#!/usr/bin/env python3
import numpy as np
import scipy.integrate as integrate
import scipy.fftpack   as fftpack
import os.path

airfoil_motions = ['M1', 'M2']

## AFRL Cases 
#afrl_grids   = ['h0', 'h1', 'h2']
#afrl_motions = ['M1', 'M2', 'M3', 'M4']
#afrl_physics = ['Re10','Re1000','ReInf']
#afrl_orders  = ['p1', 'p2', 'p3']
#afrl_times   = ['0.1', '0.01', '0.001']
#
#for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
#    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):
#
#
#        if physic == 'Re10':
#            motion1_truth = ['h2','p3','0.01']
#            motion2_truth = ['h2','p3','0.01']
#            motion3_truth = ['h2','p3','0.01']
#            motion4_truth = ['h0','p3','0.001']
#        elif physic == 'Re1000':
#            motion1_truth = ['h2','p3','0.001']
#            motion2_truth = ['h2','p3','0.001']
#            motion3_truth = ['h2','p3','0.001']
#            motion4_truth = ['h0','p2','0.01']
#        elif physic == 'ReInf':
#            motion1_truth = ['h0','p1','0.1']
#            motion2_truth = ['h0','p3','0.001']
#            motion3_truth = ['h0','p3','0.01']
#            motion4_truth = ['h0','p1','0.01']
#        else:
#            print("ERROR!!!!!")
#
#
#        if (motion == 'M1'):
#            h_ref = motion1_truth[0]
#            p_ref = motion1_truth[1]
#            t_ref = motion1_truth[2]
#        elif (motion == 'M2'):
#            h_ref = motion2_truth[0]
#            p_ref = motion2_truth[1]
#            t_ref = motion2_truth[2]
#        elif (motion == 'M3'):
#            h_ref = motion3_truth[0]
#            p_ref = motion3_truth[1]
#            t_ref = motion3_truth[2]
#        elif (motion == 'M4'):
#            h_ref = motion4_truth[0]
#            p_ref = motion4_truth[1]
#            t_ref = motion4_truth[2]
#
#        afrl_file = 'AFRL/'+motion+'-'+physic+'-'+h_ref+'-'+p_ref+'/'+'cyl-t'+t_ref
#    
#        # Load data
#        if (os.path.isfile(afrl_file)):
#            afrl_data = np.loadtxt(afrl_file, skiprows=1)
#        else:
#            afrl_data = False
#        
#        # Reported data includes initial and final time
#        #   t0 --- t1 --- t2 --- t3
#        #
#        #   e.g.
#        #   nt     = 4
#        #   nsteps = 3
#        
#            
#        if ( isinstance(afrl_data, np.ndarray) ):
#            # AFRL Data
#            afrl_Fx   = afrl_data[:,0]
#            afrl_Fy   = afrl_data[:,1]
#            afrl_Fz   = afrl_data[:,2]
#            afrl_Wint = afrl_data[:,3]
#            afrl_mass = afrl_data[:,4]
#            
#            afrl_nx = len(afrl_Fy)
#            xend    = 2.
#            afrl_dx = 2./(afrl_nx-1)
#            afrl_x  = np.linspace(0.,xend,afrl_nx)
#            
#            if (motion == 'M1'):
#                if (physic == 'ReInf'):
#                    ax_M1_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M1_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M1_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re1000'):
#                    ax_M1_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M1_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M1_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re10'):
#                    ax_M1_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M1_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M1_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#            elif (motion == 'M2'):
#                if (physic == 'ReInf'):
#                    ax_M2_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M2_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M2_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re1000'):
#                    ax_M2_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M2_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M2_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re10'):
#                    ax_M2_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M2_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M2_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#            elif (motion == 'M3'):
#                if (physic == 'ReInf'):
#                    ax_M3_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M3_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M3_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re1000'):
#                    ax_M3_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M3_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M3_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re10'):
#                    ax_M3_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M3_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M3_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#            elif (motion == 'M4'):
#                if (physic == 'ReInf'):
#                    ax_M4_1.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M4_2.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M4_3.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re1000'):
#                    ax_M4_4.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M4_5.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M4_6.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')
#                elif (physic == 'Re10'):
#                    ax_M4_7.plot(afrl_x,afrl_Fx,  'b',linewidth=1.0,label='AFRL')
#                    ax_M4_8.plot(afrl_x,afrl_Fy,  'b',linewidth=1.0,label='AFRL')
#                    ax_M4_9.plot(afrl_x,afrl_Wint,'b',linewidth=1.0,label='AFRL')





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
        ku_file = 'KU/p2_translational.txt'

    elif motion == "M2":
        ku_file = 'KU/p2_translationalPlusRotational.txt'


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

        # University of Kansas Data
        ku_t    = ku_data[:,0]
        ku_Fx   = ku_data[:,1]
        ku_Fy   = ku_data[:,2]
        ku_Wint = ku_data[:,3]
        
        KU_xI[imotion] = integrate.simps(ku_Fx,  x=ku_t)
        KU_yI[imotion] = integrate.simps(ku_Fy,  x=ku_t)
        KU_W[ imotion] = integrate.simps(ku_Wint,x=ku_t)







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



print(" ")
print("Org.  ", "X-Impulse", "Y-Impulse", "Work")
for imotion,motion in zip(range(len(airfoil_motions)),airfoil_motions):
    print("---------------------------------------------------------------- ")
    print(" ")
    print(motion)
    #print("AFRL:  ", afrl_xI[imotion], afrl_yI[imotion], afrl_W[imotion])
    print("UMich: ",   UM_xI[imotion],   UM_yI[imotion],   UM_W[imotion])
    print("KU:    ",   KU_xI[imotion],   KU_yI[imotion],   KU_W[imotion])




