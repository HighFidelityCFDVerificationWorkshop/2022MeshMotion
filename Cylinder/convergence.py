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



# UM Cases 
case_grids          = ['h0', 'h1', 'h2']
case_grid_sizes     = [20, 80, 320]
case_motions        = ['M1', 'M2', 'M3', 'M4']
case_orders         = ['p1','p2','p3','p4']
# Re3=Euler, Re2=Re1000, Re1=Re10
case_physics        = ['Re3','Re2','Re1']
case_order_integers = [1, 2, 3, 4]
case_times          = ['cyl0', 'cyl1', 'cyl2', 'cyl3']

# Storage (motions,physics,h,p,t)
xI = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders),len(case_times)))
yI = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders),len(case_times)))
zI = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders),len(case_times)))
W  = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders),len(case_times)))
m  = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders),len(case_times)))

for imotion,motion in zip(range(len(case_motions)),case_motions):
    for iphysic,physic in zip(range(len(case_physics)),case_physics):
        for igrid,grid in zip(range(len(case_grids)),case_grids):
            for iorder,order in zip(range(len(case_orders)),case_orders):
                for itime,time in zip(range(len(case_times)),case_times):

                    case_file = 'UMnew/Quad-'+motion+'-'+physic+'-'+grid+'-'+order+'/'+time+'_TimeHist.txt'

                    # Load data
                    if (os.path.isfile(case_file)):
                        case_data = np.loadtxt(case_file, skiprows=37)
                    else:
                        print('File not found:' + case_file)
                        case_data = False
                    
                    # Reported data includes initial and final time
                    #   t0 --- t1 --- t2 --- t3
                    #
                    #   e.g.
                    #   nt     = 4
                    #   nsteps = 3
                        
                    if ( isinstance(case_data, np.ndarray) ):
                        case_t    = case_data[:,1]
                        case_Fx   = case_data[:,4]
                        case_Fy   = case_data[:,5]
                        case_Wint = case_data[:,7]

                        xI[imotion,iphysic,igrid,iorder,itime] = integrate.simps(case_Fx,  x=case_t)
                        yI[imotion,iphysic,igrid,iorder,itime] = integrate.simps(case_Fy,  x=case_t)
                        W[ imotion,iphysic,igrid,iorder,itime] = integrate.simps(case_Wint,x=case_t)




itime = 3
order = 'p4'
iorder = 3
for imotion,motion in zip(range(len(case_motions)),case_motions):
    for iphysic,physic in zip(range(len(case_physics)),case_physics):

        sndofs = np.zeros(len(case_grids))
        for igrid,grid in zip(range(len(case_grids)),case_grids):
            # QUAD
            sndofs[igrid] = case_grid_sizes[igrid]*(case_order_integers[iorder]+1)**2.
            
        color = 'yo-'
    
        if (motion == 'M1'):
            if (physic == 'Re3'):
                ax_M1_1.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M1_2.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M1_3.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
            elif (physic == 'Re2'):
                ax_M1_4.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M1_5.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M1_6.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
            elif (physic == 'Re1'):
                ax_M1_7.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M1_8.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M1_9.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
        elif (motion == 'M2'):
            if (physic == 'Re3'):
                ax_M2_1.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M2_2.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M2_3.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
            elif (physic == 'Re2'):
                ax_M2_4.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M2_5.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M2_6.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')

                # New: M2, Re1000
                ax_M2_6.semilogx(sndofs[1],-1.124887156483727E-01,'go-',linewidth=1.0,label='UM NEW')

            elif (physic == 'Re1'):
                ax_M2_7.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M2_8.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M2_9.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
        elif (motion == 'M3'):
            if (physic == 'Re3'):
                ax_M3_1.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M3_2.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M3_3.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
            elif (physic == 'Re2'):
                ax_M3_4.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M3_5.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M3_6.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
            elif (physic == 'Re1'):
                ax_M3_7.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M3_8.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M3_9.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
        elif (motion == 'M4'):
            if (physic == 'Re3'):
                ax_M4_1.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M4_2.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M4_3.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
            elif (physic == 'Re2'):
                ax_M4_4.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M4_5.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M4_6.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
            elif (physic == 'Re1'):
                ax_M4_7.semilogx(sndofs,xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M4_8.semilogx(sndofs,yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')
                ax_M4_9.semilogx(sndofs,W[ imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='UM P4')


# AFRL Cases 
afrl_grids          = ['h0', 'h1', 'h2']
afrl_grid_sizes     = [20, 80, 320]
afrl_motions        = ['M1', 'M2', 'M3', 'M4']
afrl_physics        = ['ReInf','Re1000','Re10']
afrl_orders         = ['p1','p2','p3']
afrl_order_integers = [1, 2, 3, 4]
afrl_times          = ['0.1','0.01','0.001']

# Storage (motions,physics,h,p,t)
afrl_xI = np.zeros((len(afrl_motions),len(afrl_physics),len(afrl_grids),len(afrl_orders),len(afrl_times)))
afrl_yI = np.zeros((len(afrl_motions),len(afrl_physics),len(afrl_grids),len(afrl_orders),len(afrl_times)))
afrl_zI = np.zeros((len(afrl_motions),len(afrl_physics),len(afrl_grids),len(afrl_orders),len(afrl_times)))
afrl_W  = np.zeros((len(afrl_motions),len(afrl_physics),len(afrl_grids),len(afrl_orders),len(afrl_times)))
afrl_m  = np.zeros((len(afrl_motions),len(afrl_physics),len(afrl_grids),len(afrl_orders),len(afrl_times)))


for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):
        for igrid,grid in zip(range(len(afrl_grids)),afrl_grids):
            for iorder,order in zip(range(len(afrl_orders)),afrl_orders):
                for itime,time in zip(range(len(afrl_times)),afrl_times):


                    afrl_file = 'AFRL/'+motion+'-'+physic+'-'+grid+'-'+order+'/cyl-t'+time
                
                    # Load data
                    if (os.path.isfile(afrl_file)):
                        afrl_data = np.loadtxt(afrl_file, skiprows=1, ndmin=2)
                    else:
                        print('File not found: ' + afrl_file)
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

                        if afrl_nx == 21 or afrl_nx == 201 or afrl_nx == 2001:

                            afrl_dx = 2./(afrl_nx-1)
                            afrl_x  = np.linspace(0.,xend,afrl_nx)

                            afrl_xI[imotion,iphysic,igrid,iorder,itime] = integrate.simps(afrl_Fx,  dx=afrl_dx)
                            afrl_yI[imotion,iphysic,igrid,iorder,itime] = integrate.simps(afrl_Fy,  dx=afrl_dx)
                            afrl_zI[imotion,iphysic,igrid,iorder,itime] = integrate.simps(afrl_Fz,  dx=afrl_dx)
                            afrl_W[ imotion,iphysic,igrid,iorder,itime] = integrate.simps(afrl_Wint,dx=afrl_dx)
                            afrl_m[ imotion,iphysic,igrid,iorder,itime] = integrate.simps(afrl_mass-afrl_mass[0],dx=afrl_dx)
                        else:
                            print("bad data from case: " + afrl_file)




#itime = len(afrl_times)-1
itime = 2
iorder = 2
order = 'p3'
for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):

        sndofs = np.zeros(len(afrl_grids))
        for igrid,grid in zip(range(len(afrl_grids)),afrl_grids):
            sndofs[igrid] = afrl_grid_sizes[igrid]*(afrl_order_integers[iorder]+1)**2.
            
        color = 'bo-'
    
        if (motion == 'M1'):
            if (physic == 'ReInf'):
                ax_M1_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M1_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M1_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re1000'):
                ax_M1_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M1_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M1_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re10'):
                ax_M1_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M1_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M1_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')
        elif (motion == 'M2'):
            if (physic == 'ReInf'):
                ax_M2_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M2_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M2_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re1000'):
                ax_M2_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M2_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M2_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re10'):
                ax_M2_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M2_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M2_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')
        elif (motion == 'M3'):
            if (physic == 'ReInf'):
                ax_M3_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M3_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M3_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re1000'):
                ax_M3_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M3_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M3_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re10'):
                ax_M3_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M3_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M3_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')
        elif (motion == 'M4'):
            if (physic == 'ReInf'):
                ax_M4_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M4_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M4_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re1000'):
                ax_M4_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M4_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M4_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')
            elif (physic == 'Re10'):
                ax_M4_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M4_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label='AFRL P3')
                ax_M4_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label='AFRL P3')








# UCB Cases 
case_grids          = ['ref0', 'ref1', 'ref2', 'ref3']
case_grid_sizes     = [20, 80, 320, 1280]
case_motions        = ['M1', 'M2', 'M3', 'M4']
case_orders         = ['p1', 'p2', 'p3', 'p4']
case_physics        = ['ReInf','Re1000','Re10']
case_order_integers = [1, 2, 3, 4]

# Storage (motions,physics,h,p,t)
xI = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders)))
yI = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders)))
zI = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders)))
W  = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders)))
m  = np.zeros((len(case_motions),len(case_physics),len(case_grids),len(case_orders)))

for imotion,motion in zip(range(len(case_motions)),case_motions):
    for iphysic,physic in zip(range(len(case_physics)),case_physics):
        for igrid,grid in zip(range(len(case_grids)),case_grids):
            for iorder,order in zip(range(len(case_orders)),case_orders):

                case_file = 'UCB/'+motion+physic+'_'+grid+order+'.dat'

                # Load data
                if (os.path.isfile(case_file)):
                    case_data = np.loadtxt(case_file, skiprows=1)
                else:
                    print(case_file)
                    print('UCB data not found!')
                    case_data = False


                # Reported data includes initial and final time
                #   t0 --- t1 --- t2 --- t3
                #
                #   e.g.
                #   nt     = 4
                #   nsteps = 3
                    
                if ( isinstance(case_data, np.ndarray) ):
                    # U.C. Berkeley Data
                    case_t    = case_data[:,0]
                    case_Fx   = case_data[:,1]
                    case_Fy   = case_data[:,2]
                    case_Wint = case_data[:,3]



                    xI[imotion,iphysic,igrid,iorder] = integrate.simps(case_Fx,  x=case_t)
                    yI[imotion,iphysic,igrid,iorder] = integrate.simps(case_Fy,  x=case_t)
                    W[ imotion,iphysic,igrid,iorder] = integrate.simps(case_Wint,x=case_t)



iorder = 2
order = 'p3'
for imotion,motion in zip(range(len(case_motions)),case_motions):
    for iphysic,physic in zip(range(len(case_physics)),case_physics):


        #sndofs = np.zeros(len(case_grids))
        #for igrid,grid in zip(range(len(case_grids)),case_grids):
        #    # Tri
        #    sndofs[igrid] = case_grid_sizes[igrid]*(case_order_integers[iorder]+1)**2.

        #sndofs = np.zeros(len(case_grids))
        #for igrid,grid in zip(range(len(case_grids)),case_grids):
        #    sndofs[igrid] = case_grid_sizes[igrid]*(afrl_order_integers[iorder]+1)**2.

        ##sndofs = [1., 2., 3., 4.]

        sndofs = [40., 160., 640., 2560.]
            
        color = 'ro-'
    
        if (motion == 'M1'):
            if (physic == 'ReInf'):
                ax_M1_1.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M1_2.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M1_3.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re1000'):
                ax_M1_4.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M1_5.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M1_6.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re10'):
                ax_M1_7.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M1_8.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M1_9.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
        elif (motion == 'M2'):
            if (physic == 'ReInf'):
                ax_M2_1.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M2_2.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M2_3.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re1000'):
                ax_M2_4.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M2_5.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M2_6.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re10'):
                if order == 'p4':
                    ax_M2_7.plot(sndofs[:-1],xI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M2_8.plot(sndofs[:-1],yI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M2_9.plot(sndofs[:-1],W[ imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                else:
                    ax_M2_7.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M2_8.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M2_9.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
        elif (motion == 'M3'):
            if (physic == 'ReInf'):
                ax_M3_1.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M3_2.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M3_3.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re1000'):
                if order == 'p4':
                    ax_M3_4.plot(sndofs[:-1],xI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M3_5.plot(sndofs[:-1],yI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M3_6.plot(sndofs[:-1],W[ imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                else:
                    ax_M3_4.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M3_5.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M3_6.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re10'):
                if order == 'p1' or order == 'p2' or order == 'p4':
                    ax_M3_7.plot(sndofs[:-1],xI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M3_8.plot(sndofs[:-1],yI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M3_9.plot(sndofs[:-1],W[ imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                else:
                    ax_M3_7.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M3_8.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M3_9.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
        elif (motion == 'M4'):
            if (physic == 'ReInf'):
                ax_M4_1.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M4_2.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                ax_M4_3.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re1000'):
                if order == 'p4':
                    ax_M4_4.plot(sndofs[:-1],xI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M4_5.plot(sndofs[:-1],yI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M4_6.plot(sndofs[:-1],W[ imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                else:
                    ax_M4_4.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M4_5.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M4_6.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
            elif (physic == 'Re10'):
                if order == 'p4':
                    ax_M4_7.plot(sndofs[:-2],xI[imotion,iphysic,:-2,iorder],color,linewidth=1.0,label=order)
                    ax_M4_8.plot(sndofs[:-2],yI[imotion,iphysic,:-2,iorder],color,linewidth=1.0,label=order)
                    ax_M4_9.plot(sndofs[:-2],W[ imotion,iphysic,:-2,iorder],color,linewidth=1.0,label=order)
                elif order == 'p1' or order == 'p2':
                    ax_M4_7.plot(sndofs[:-1],xI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M4_8.plot(sndofs[:-1],yI[imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                    ax_M4_9.plot(sndofs[:-1],W[ imotion,iphysic,:-1,iorder],color,linewidth=1.0,label=order)
                else:
                    ax_M4_7.plot(sndofs,xI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M4_8.plot(sndofs,yI[imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')
                    ax_M4_9.plot(sndofs,W[ imotion,iphysic,:,iorder],color,linewidth=1.0,label='UCB P3')























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

ax_M1_1.set_xlabel("Spatial DOF's")
ax_M1_2.set_xlabel("Spatial DOF's")
ax_M1_3.set_xlabel("Spatial DOF's")
ax_M1_4.set_xlabel("Spatial DOF's")
ax_M1_5.set_xlabel("Spatial DOF's")
ax_M1_6.set_xlabel("Spatial DOF's")
ax_M1_7.set_xlabel("Spatial DOF's")
ax_M1_8.set_xlabel("Spatial DOF's")
ax_M1_9.set_xlabel("Spatial DOF's")

ax_M2_1.set_xlabel("Spatial DOF's")
ax_M2_2.set_xlabel("Spatial DOF's")
ax_M2_3.set_xlabel("Spatial DOF's")
ax_M2_4.set_xlabel("Spatial DOF's")
ax_M2_5.set_xlabel("Spatial DOF's")
ax_M2_6.set_xlabel("Spatial DOF's")
ax_M2_7.set_xlabel("Spatial DOF's")
ax_M2_8.set_xlabel("Spatial DOF's")
ax_M2_9.set_xlabel("Spatial DOF's")

ax_M3_1.set_xlabel("Spatial DOF's")
ax_M3_2.set_xlabel("Spatial DOF's")
ax_M3_3.set_xlabel("Spatial DOF's")
ax_M3_4.set_xlabel("Spatial DOF's")
ax_M3_5.set_xlabel("Spatial DOF's")
ax_M3_6.set_xlabel("Spatial DOF's")
ax_M3_7.set_xlabel("Spatial DOF's")
ax_M3_8.set_xlabel("Spatial DOF's")
ax_M3_9.set_xlabel("Spatial DOF's")

ax_M4_1.set_xlabel("Spatial DOF's")
ax_M4_2.set_xlabel("Spatial DOF's")
ax_M4_3.set_xlabel("Spatial DOF's")
ax_M4_4.set_xlabel("Spatial DOF's")
ax_M4_5.set_xlabel("Spatial DOF's")
ax_M4_6.set_xlabel("Spatial DOF's")
ax_M4_7.set_xlabel("Spatial DOF's")
ax_M4_8.set_xlabel("Spatial DOF's")
ax_M4_9.set_xlabel("Spatial DOF's")

ax_M1_1.set_ylabel('Impulse-X (Euler)')
ax_M1_2.set_ylabel('Impulse-Y (Euler)')
ax_M1_3.set_ylabel('Work (Euler)')
ax_M1_4.set_ylabel('Impulse-X (Re = 1000)')
ax_M1_5.set_ylabel('Impulse-Y (Re = 1000)')
ax_M1_6.set_ylabel('Work (Re = 1000)')
ax_M1_7.set_ylabel('Impulse-X (Re = 10)')
ax_M1_8.set_ylabel('Impulse-Y (Re = 10)')
ax_M1_9.set_ylabel('Work (Re = 10)')

ax_M2_1.set_ylabel('Impulse-X (Euler)')
ax_M2_2.set_ylabel('Impulse-Y (Euler)')
ax_M2_3.set_ylabel('Work (Euler)')
ax_M2_4.set_ylabel('Impulse-X (Re = 1000)')
ax_M2_5.set_ylabel('Impulse-Y (Re = 1000)')
ax_M2_6.set_ylabel('Work (Re = 1000)')
ax_M2_7.set_ylabel('Impulse-X (Re = 10)')
ax_M2_8.set_ylabel('Impulse-Y (Re = 10)')
ax_M2_9.set_ylabel('Work (Re = 10)')

ax_M3_1.set_ylabel('Impulse-X (Euler)')
ax_M3_2.set_ylabel('Impulse-Y (Euler)')
ax_M3_3.set_ylabel('Work (Euler)')
ax_M3_4.set_ylabel('Impulse-X (Re = 1000)')
ax_M3_5.set_ylabel('Impulse-Y (Re = 1000)')
ax_M3_6.set_ylabel('Work (Re = 1000)')
ax_M3_7.set_ylabel('Impulse-X (Re = 10)')
ax_M3_8.set_ylabel('Impulse-Y (Re = 10)')
ax_M3_9.set_ylabel('Work (Re = 10)')

ax_M4_1.set_ylabel('Impulse-X (Euler)')
ax_M4_2.set_ylabel('Impulse-Y (Euler)')
ax_M4_3.set_ylabel('Work (Euler)')
ax_M4_4.set_ylabel('Impulse-X (Re = 1000)')
ax_M4_5.set_ylabel('Impulse-Y (Re = 1000)')
ax_M4_6.set_ylabel('Work (Re = 1000)')
ax_M4_7.set_ylabel('Impulse-X (Re = 10)')
ax_M4_8.set_ylabel('Impulse-Y (Re = 10)')
ax_M4_9.set_ylabel('Work (Re = 10)')



#ax_M1_1.set_xlim((0.,2.))
#ax_M1_2.set_xlim((0.,2.))
#ax_M1_3.set_xlim((0.,2.))
#ax_M1_4.set_xlim((0.,2.))
#ax_M1_5.set_xlim((0.,2.))
#ax_M1_6.set_xlim((0.,2.))
#ax_M1_7.set_xlim((0.,2.))
#ax_M1_8.set_xlim((0.,2.))
#ax_M1_9.set_xlim((0.,2.))
##
#ax_M2_1.set_xlim((0.,2.))
#ax_M2_2.set_xlim((0.,2.))
#ax_M2_3.set_xlim((0.,2.))
#ax_M2_4.set_xlim((0.,2.))
#ax_M2_5.set_xlim((0.,2.))
#ax_M2_6.set_xlim((0.,2.))
#ax_M2_7.set_xlim((0.,2.))
#ax_M2_8.set_xlim((0.,2.))
#ax_M2_9.set_xlim((0.,2.))
##
#ax_M3_1.set_xlim((0.,2.))
#ax_M3_2.set_xlim((0.,2.))
#ax_M3_3.set_xlim((0.,2.))
#ax_M3_4.set_xlim((0.,2.))
#ax_M3_5.set_xlim((0.,2.))
#ax_M3_6.set_xlim((0.,2.))
#ax_M3_7.set_xlim((0.,2.))
#ax_M3_8.set_xlim((0.,2.))
#ax_M3_9.set_xlim((0.,2.))
##
#ax_M4_1.set_xlim((0.,2.))
#ax_M4_2.set_xlim((0.,2.))
#ax_M4_3.set_xlim((0.,2.))
#ax_M4_4.set_xlim((0.,2.))
#ax_M4_5.set_xlim((0.,2.))
#ax_M4_6.set_xlim((0.,2.))
#ax_M4_7.set_xlim((0.,2.))
#ax_M4_8.set_xlim((0.,2.))
#ax_M4_9.set_xlim((0.,2.))

ax_M1_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M1_4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M1_7.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M2_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M2_2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M2_4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M2_5.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M2_7.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M2_8.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M3_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M3_2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M3_4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M3_5.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M3_7.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_M3_8.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


ax_M1_1.set_ylim((-0.1e-08,0.1e-08))
ax_M1_2.set_ylim((-0.0118,-0.0116))
ax_M1_3.set_ylim((-0.0088,-0.0086))
ax_M1_4.set_ylim((-0.1e-08,0.1e-08))
ax_M1_5.set_ylim((-0.0125,-0.0120))
ax_M1_6.set_ylim((-0.0088,-0.008))
ax_M1_7.set_ylim((-0.1e-08,0.1e-08))
ax_M1_8.set_ylim((-0.014,-0.012))
ax_M1_9.set_ylim((-0.01025,-0.010))

ax_M2_1.set_ylim((-0.1e-8,0.1e-8))
ax_M2_2.set_ylim((-0.1e-8,0.1e-8))
ax_M2_3.set_ylim((-0.1e-6,0.1e-6))
ax_M2_4.set_ylim((-0.1e-8,0.1e-8))
ax_M2_5.set_ylim((-0.1e-8,0.1e-8))
ax_M2_6.set_ylim((-0.12,-0.1))
ax_M2_7.set_ylim((-0.1e-8,0.1e-8))
ax_M2_8.set_ylim((-0.1e-8,0.1e-8))
ax_M2_9.set_ylim((-0.178,-0.174))

ax_M3_1.set_ylim((-0.1e-08,0.1e-08))
ax_M3_2.set_ylim((-0.1e-08,0.1e-08))
ax_M3_3.set_ylim((-2.e-4,-1.e-5))
ax_M3_4.set_ylim((-0.1e-08,0.1e-08))
ax_M3_5.set_ylim((-0.1e-08,0.1e-08))
ax_M3_6.set_ylim((-0.00044,-0.000425))
ax_M3_7.set_ylim((-0.1e-08,0.1e-08))
ax_M3_8.set_ylim((-0.1e-08,0.1e-08))
ax_M3_9.set_ylim((-0.035,-0.0345))

ax_M4_1.set_ylim((-0.006,-0.005))
ax_M4_2.set_ylim((-0.010,-0.008))
ax_M4_3.set_ylim((-0.006,-0.004))
ax_M4_4.set_ylim((-0.006,-0.004))
ax_M4_5.set_ylim((-0.009,-0.006))
ax_M4_6.set_ylim((-0.1284,-0.128))
ax_M4_7.set_ylim((-0.008,-0.004))
ax_M4_8.set_ylim((-0.0065,-0.006))
ax_M4_9.set_ylim((-0.17,-0.165))






M1.savefig('Cylinder_convergence_motion1.png', bbox_inches='tight', dpi=800)
M2.savefig('Cylinder_convergence_motion2.png', bbox_inches='tight', dpi=800)
M3.savefig('Cylinder_convergence_motion3.png', bbox_inches='tight', dpi=800)
M4.savefig('Cylinder_convergence_motion4.png', bbox_inches='tight', dpi=800)

plt.show()




