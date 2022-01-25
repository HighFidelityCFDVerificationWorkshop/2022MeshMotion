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



M5 = plt.figure(figsize=(10,9))
M5.set_facecolor('White')
ax_M5_1 = M5.add_subplot(331)
ax_M5_2 = M5.add_subplot(332)
ax_M5_3 = M5.add_subplot(333)
ax_M5_4 = M5.add_subplot(334)
ax_M5_5 = M5.add_subplot(335)
ax_M5_6 = M5.add_subplot(336)
ax_M5_7 = M5.add_subplot(337)
ax_M5_8 = M5.add_subplot(338)
ax_M5_9 = M5.add_subplot(339)

M6 = plt.figure(figsize=(10,9))
M6.set_facecolor('White')
ax_M6_1 = M6.add_subplot(331)
ax_M6_2 = M6.add_subplot(332)
ax_M6_3 = M6.add_subplot(333)
ax_M6_4 = M6.add_subplot(334)
ax_M6_5 = M6.add_subplot(335)
ax_M6_6 = M6.add_subplot(336)
ax_M6_7 = M6.add_subplot(337)
ax_M6_8 = M6.add_subplot(338)
ax_M6_9 = M6.add_subplot(339)

M7 = plt.figure(figsize=(10,9))
M7.set_facecolor('White')
ax_M7_1 = M7.add_subplot(331)
ax_M7_2 = M7.add_subplot(332)
ax_M7_3 = M7.add_subplot(333)
ax_M7_4 = M7.add_subplot(334)
ax_M7_5 = M7.add_subplot(335)
ax_M7_6 = M7.add_subplot(336)
ax_M7_7 = M7.add_subplot(337)
ax_M7_8 = M7.add_subplot(338)
ax_M7_9 = M7.add_subplot(339)

M8 = plt.figure(figsize=(10,9))
M8.set_facecolor('White')
ax_M8_1 = M8.add_subplot(331)
ax_M8_2 = M8.add_subplot(332)
ax_M8_3 = M8.add_subplot(333)
ax_M8_4 = M8.add_subplot(334)
ax_M8_5 = M8.add_subplot(335)
ax_M8_6 = M8.add_subplot(336)
ax_M8_7 = M8.add_subplot(337)
ax_M8_8 = M8.add_subplot(338)
ax_M8_9 = M8.add_subplot(339)





# AFRL Cases 
afrl_grids          = ['h0', 'h1', 'h2']
afrl_grid_sizes     = [20, 80, 320]
afrl_motions        = ['M1', 'M2', 'M3', 'M4']
afrl_physics        = ['ReInf','Re1000','Re10']
afrl_orders         = ['p1', 'p2', 'p3']
afrl_order_integers = [1, 2, 3]
afrl_times          = ['0.1', '0.01', '0.001']

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


                    afrl_file = motion+'-'+physic+'-'+grid+'-'+order+'/cyl-t'+time
                
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
for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):
        for iorder,order in zip(range(len(afrl_orders)),afrl_orders):


            sndofs = np.zeros(len(afrl_grids))
            for igrid,grid in zip(range(len(afrl_grids)),afrl_grids):
                sndofs[igrid] = afrl_grid_sizes[igrid]*(afrl_order_integers[iorder]+1)**2.
                
            if order == 'p1':
                color = 'ro-'
            elif order == 'p2':
                color = 'bo-'
            elif order == 'p3':
                color = 'go-'
            elif order == 'p4':
                color = 'ko-'

        
            if (motion == 'M1'):
                if (physic == 'ReInf'):
                    ax_M1_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M1_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M1_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    if order == 'p1':
                        ax_M1_4.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                        ax_M1_5.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                        ax_M1_6.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label=order)
                    else:
                        ax_M1_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                        ax_M1_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                        ax_M1_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M1_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M1_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M1_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)
            elif (motion == 'M2'):
                if (physic == 'ReInf'):
                    ax_M2_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M2_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M2_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    ax_M2_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M2_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M2_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M2_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M2_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M2_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)
            elif (motion == 'M3'):
                if (physic == 'ReInf'):
                    ax_M3_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M3_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M3_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    ax_M3_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M3_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M3_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M3_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M3_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M3_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)
            elif (motion == 'M4'):
                if (physic == 'ReInf'):
                    ax_M4_1.semilogx(sndofs[:-1],afrl_xI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M4_2.semilogx(sndofs[:-1],afrl_yI[imotion,iphysic,:-1,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M4_3.semilogx(sndofs[:-1],afrl_W[ imotion,iphysic,:-1,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    ax_M4_4.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M4_5.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M4_6.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M4_7.semilogx(sndofs,afrl_xI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M4_8.semilogx(sndofs,afrl_yI[imotion,iphysic,:,iorder,itime],color,linewidth=1.0,label=order)
                    ax_M4_9.semilogx(sndofs,afrl_W[ imotion,iphysic,:,iorder,itime], color,linewidth=1.0,label=order)




# Truth data storage (motions,physics)
truth_xI = np.zeros((len(afrl_motions),len(afrl_physics)))
truth_yI = np.zeros((len(afrl_motions),len(afrl_physics)))
truth_zI = np.zeros((len(afrl_motions),len(afrl_physics)))
truth_W  = np.zeros((len(afrl_motions),len(afrl_physics)))

for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):

        if motion == "ReInf":
            igrid = 1
            iorder = 2
            itime = 2
            truth_xI[imotion,iphysic] = afrl_xI[imotion,iphysic,igrid,iorder,itime]
            truth_yI[imotion,iphysic] = afrl_yI[imotion,iphysic,igrid,iorder,itime]
            truth_zI[imotion,iphysic] = afrl_zI[imotion,iphysic,igrid,iorder,itime]
            truth_W[ imotion,iphysic] = afrl_W[ imotion,iphysic,igrid,iorder,itime]
        else:
            igrid = 2
            iorder = 2
            itime = 2
            truth_xI[imotion,iphysic] = afrl_xI[imotion,iphysic,igrid,iorder,itime]
            truth_yI[imotion,iphysic] = afrl_yI[imotion,iphysic,igrid,iorder,itime]
            truth_zI[imotion,iphysic] = afrl_zI[imotion,iphysic,igrid,iorder,itime]
            truth_W[ imotion,iphysic] = afrl_W[ imotion,iphysic,igrid,iorder,itime]
        




itime = 2
for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):
        for iorder,order in zip(range(len(afrl_orders)),afrl_orders):


            sndofs = np.zeros(len(afrl_grids))
            for igrid,grid in zip(range(len(afrl_grids)),afrl_grids):
                sndofs[igrid] = afrl_grid_sizes[igrid]*(afrl_order_integers[iorder]+1)**2.
                

            if order == 'p1':
                color = 'ro-'
            elif order == 'p2':
                color = 'bo-'
            elif order == 'p3':
                color = 'go-'
            elif order == 'p4':
                color = 'ko-'

        
            if (motion == 'M1'):
                if (physic == 'ReInf'):
                    ax_M5_1.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M5_2.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M5_3.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    if order == 'p1':
                        ax_M5_4.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                        ax_M5_5.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                        ax_M5_6.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    else:
                        ax_M5_4.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                        ax_M5_5.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                        ax_M5_6.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M5_7.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M5_8.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M5_9.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
            elif (motion == 'M2'):
                if (physic == 'ReInf'):
                    ax_M6_1.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M6_2.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M6_3.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    ax_M6_4.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M6_5.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M6_6.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M6_7.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M6_8.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M6_9.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
            elif (motion == 'M3'):
                if (physic == 'ReInf'):
                    ax_M7_1.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M7_2.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M7_3.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    ax_M7_4.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M7_5.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M7_6.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M7_7.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M7_8.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M7_9.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
            elif (motion == 'M4'):
                if (physic == 'ReInf'):
                    ax_M8_1.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M8_2.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M8_3.loglog(1./np.sqrt(sndofs[:-1]),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:-1,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re1000'):
                    ax_M8_4.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M8_5.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M8_6.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                elif (physic == 'Re10'):
                    ax_M8_7.loglog(1./np.sqrt(sndofs),np.abs(truth_xI[imotion,iphysic] - afrl_xI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M8_8.loglog(1./np.sqrt(sndofs),np.abs(truth_yI[imotion,iphysic] - afrl_yI[imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)
                    ax_M8_9.loglog(1./np.sqrt(sndofs),np.abs(truth_W[ imotion,iphysic] - afrl_W[ imotion,iphysic,:,iorder,itime]),color,linewidth=1.0,label=order)


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

ax_M1_1.set_xlabel("Shatial DOF's")
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
ax_M1_2.set_ylim((-0.012,-0.011))
ax_M1_3.set_ylim((-0.009,-0.008))
ax_M1_4.set_ylim((-0.1e-08,0.1e-08))
ax_M1_5.set_ylim((-0.014,-0.01))
ax_M1_6.set_ylim((-0.009,-0.008))
ax_M1_7.set_ylim((-0.1e-08,0.1e-08))
ax_M1_8.set_ylim((-0.015,-0.005))
ax_M1_9.set_ylim((-0.012,-0.005))

ax_M2_1.set_ylim((-0.1e-8,0.1e-8))
ax_M2_2.set_ylim((-0.1e-8,0.1e-8))
ax_M2_3.set_ylim((-0.1e-8,0.1e-8))
ax_M2_4.set_ylim((-0.1e-8,0.1e-8))
ax_M2_5.set_ylim((-0.1e-8,0.1e-8))
ax_M2_6.set_ylim((-0.12,-0.1))
ax_M2_7.set_ylim((-0.1e-8,0.1e-8))
ax_M2_8.set_ylim((-0.1e-8,0.1e-8))
ax_M2_9.set_ylim((-0.20,-0.15))

ax_M3_1.set_ylim((-0.1e-08,0.1e-08))
ax_M3_2.set_ylim((-0.1e-08,0.1e-08))
ax_M3_3.set_ylim((-1.e-5,-2.e-4))
ax_M3_4.set_ylim((-0.1e-08,0.1e-08))
ax_M3_5.set_ylim((-0.1e-08,0.1e-08))
ax_M3_6.set_ylim((-0.0006,-0.0003))
ax_M3_7.set_ylim((-0.1e-08,0.1e-08))
ax_M3_8.set_ylim((-0.1e-08,0.1e-08))
ax_M3_9.set_ylim((-0.04,-0.03))

ax_M4_1.set_ylim((-0.006,-0.005))
ax_M4_2.set_ylim((-0.010,-0.008))
ax_M4_3.set_ylim((-0.006,-0.004))
ax_M4_4.set_ylim((-0.006,-0.004))
ax_M4_5.set_ylim((-0.009,-0.006))
ax_M4_6.set_ylim((-0.130,-0.12))
ax_M4_7.set_ylim((-0.008,-0.004))
ax_M4_8.set_ylim((-0.008,-0.004))
ax_M4_9.set_ylim((-0.20,-0.12))



M1.savefig('AFRL_convergence_motion1.png', bbox_inches='tight', dpi=800)
M2.savefig('AFRL_convergence_motion2.png', bbox_inches='tight', dpi=800)
M3.savefig('AFRL_convergence_motion3.png', bbox_inches='tight', dpi=800)
M4.savefig('AFRL_convergence_motion4.png', bbox_inches='tight', dpi=800)

plt.show()




