#!/usr/bin/env python3
import numpy as np
import scipy.integrate as integrate
import scipy.fftpack   as fftpack
import os.path

# AFRL Cases 
afrl_motions = ['M1', 'M2', 'M3', 'M4']
afrl_physics = ['ReInf','Re1000','Re10']

# Storage (motions,physics)
afrl_xI = np.zeros((4,3))
afrl_yI = np.zeros((4,3))
afrl_zI = np.zeros((4,3))
afrl_W  = np.zeros((4,3))
afrl_m  = np.zeros((4,3))

for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):

        if physic == 'Re10':
            motion1_truth = ['h2','p3','0.01']
            motion2_truth = ['h2','p3','0.01']
            motion3_truth = ['h2','p3','0.01']
            motion4_truth = ['h2','p3','0.01']
        elif physic == 'Re1000':
            motion1_truth = ['h2','p3','0.01']
            motion2_truth = ['h2','p3','0.01']
            motion3_truth = ['h2','p3','0.01']
            motion4_truth = ['h2','p3','0.01']
        elif physic == 'ReInf':
            motion1_truth = ['h1','p3','0.01']
            motion2_truth = ['h2','p3','0.01']
            motion3_truth = ['h2','p3','0.01']
            motion4_truth = ['h1','p3','0.001']
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
            afrl_dx = 2./float(afrl_nx-1)
            afrl_x  = np.linspace(0.,xend,afrl_nx)
            
            afrl_xI[imotion,iphysic] = integrate.simps(afrl_Fx,  dx=afrl_dx)
            afrl_yI[imotion,iphysic] = integrate.simps(afrl_Fy,  dx=afrl_dx)
            afrl_zI[imotion,iphysic] = integrate.simps(afrl_Fz,  dx=afrl_dx)
            afrl_W[ imotion,iphysic] = integrate.simps(afrl_Wint,dx=afrl_dx)
            afrl_m[ imotion,iphysic] = integrate.simps(afrl_mass-afrl_mass[0],dx=afrl_dx)





# UM Cases 
UM_motions = ['M1', 'M2', 'M3', 'M4']
UM_physics = ['Re3','Re2','Re1']

# Storage (motions,physics)
UM_xI = np.zeros((4,3))
UM_yI = np.zeros((4,3))
UM_zI = np.zeros((4,3))
UM_W  = np.zeros((4,3))
UM_m  = np.zeros((4,3))

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
            um_t    = um_data[:,1]
            um_Fx   = um_data[:,4]
            um_Fy   = um_data[:,5]
            um_Wint = um_data[:,7]
            
            UM_xI[imotion,iphysic] = integrate.simps(um_Fx,  x=um_t)
            UM_yI[imotion,iphysic] = integrate.simps(um_Fy,  x=um_t)
            UM_W[ imotion,iphysic] = integrate.simps(um_Wint,x=um_t)
        


# KU Cases 
KU_motions = ['M1', 'M2', 'M3', 'M4']
KU_physics = ['ReInf','Re1000','Re10']

# Storage (motions,physics)
KU_xI = np.zeros((4,3))
KU_yI = np.zeros((4,3))
KU_zI = np.zeros((4,3))
KU_W  = np.zeros((4,3))
KU_m  = np.zeros((4,3))

for imotion,motion in zip(range(len(KU_motions)),KU_motions):
    for iphysic,physic in zip(range(len(KU_physics)),KU_physics):

        if motion == "M1":
            if physic == 'Re10':
                ku_file = 'KU/p3_translational_implicit.txt'
            elif physic == 'ReInf':
                ku_file = 'KU/p3_translational_Euler_implicit.txt'
            else:
                ku_file = False

        elif motion == "M2":
            if physic == 'Re10':
                ku_file = 'KU/p3_rotational_implicit.txt'
            elif physic == 'ReInf':
                ku_file = 'KU/p3_rotational_Euler_implicit.txt'
            else:
                ku_file = False

        else:
            ku_file = False


        # Load data
        if (os.path.isfile(ku_file)):
            ku_data = np.loadtxt(ku_file,delimiter=',')
        else:
            ku_data = False
    
        # Reported data does NOT include initial time. Additionally, final
        # time is duplicated.
        #   Missing --- t1 --- t2 --- ... --- t_end --- t_end
        #
        if ( isinstance(ku_data, np.ndarray) ):

            # Insert info for t=0
            ku_data = np.append([[0., 0., 0., 0.]],ku_data,axis=0)

            # Remove duplicate entry for t_end
            ku_data = np.delete(ku_data, -1, axis=0)

            # University of Michigan Data
            ku_t    = ku_data[:,0]
            ku_Fx   = ku_data[:,1]
            ku_Fy   = ku_data[:,2]
            ku_Wint = ku_data[:,3]
            
            KU_xI[imotion,iphysic] = integrate.simps(ku_Fx,  x=ku_t)
            KU_yI[imotion,iphysic] = integrate.simps(ku_Fy,  x=ku_t)
            KU_W[ imotion,iphysic] = integrate.simps(ku_Wint,x=ku_t)




print(" ")
print("Org.  ", "X-Impulse", "Y-Impulse", "Work")
for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    print("---------------------------------------------------------------- ")
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):
        print(" ")
        print(motion, physic)
        print("AFRL:  ", afrl_xI[imotion,iphysic], afrl_yI[imotion,iphysic], afrl_W[imotion,iphysic])
        print("UMich: ",   UM_xI[imotion,iphysic],   UM_yI[imotion,iphysic],   UM_W[imotion,iphysic])
        print("KU:    ",   KU_xI[imotion,iphysic],   KU_yI[imotion,iphysic],   KU_W[imotion,iphysic])



