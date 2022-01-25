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



# UCB Cases 
case_grids          = ['ref0', 'ref1', 'ref2']
case_grid_sizes     = [20, 80, 320]
case_motions        = ['M1', 'M2']
case_orders         = ['p1', 'p2', 'p3']
case_order_integers = [1, 2, 3, 4]

# Storage (motions,physics,h,p,t)
xI = np.zeros((len(case_motions),len(case_grids),len(case_orders)))
yI = np.zeros((len(case_motions),len(case_grids),len(case_orders)))
zI = np.zeros((len(case_motions),len(case_grids),len(case_orders)))
W  = np.zeros((len(case_motions),len(case_grids),len(case_orders)))
m  = np.zeros((len(case_motions),len(case_grids),len(case_orders)))

for imotion,motion in zip(range(len(case_motions)),case_motions):
    for igrid,grid in zip(range(len(case_grids)),case_grids):
        for iorder,order in zip(range(len(case_orders)),case_orders):

            case_file = motion+'_'+grid+order+'.dat'

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

                xI[imotion,igrid,iorder] = integrate.simps(case_Fx,  x=case_t)
                yI[imotion,igrid,iorder] = integrate.simps(case_Fy,  x=case_t)
                W[ imotion,igrid,iorder] = integrate.simps(case_Wint,x=case_t)




for imotion,motion in zip(range(len(case_motions)),case_motions):
    for iorder,order in zip(range(len(case_orders)),case_orders):


        #sndofs = np.zeros(len(case_grids))
        #for igrid,grid in zip(range(len(case_grids)),case_grids):
        #    # Tri
        #    sndofs[igrid] = case_grid_sizes[igrid]*(case_order_integers[iorder]+1)**2.

        sndofs = [1., 2., 3.]
            
        if order == 'p1':
            color = 'ro-'
        elif order == 'p2':
            color = 'bo-'
        elif order == 'p3':
            color = 'go-'
        elif order == 'p4':
            color = 'ko-'
    
        if (motion == 'M1'):
                ax_M1_1.plot(sndofs,xI[imotion,:,iorder],color,linewidth=1.0,label=order)
                ax_M1_2.plot(sndofs,yI[imotion,:,iorder],color,linewidth=1.0,label=order)
                ax_M1_3.plot(sndofs,W[ imotion,:,iorder],color,linewidth=1.0,label=order)
        elif (motion == 'M2'):
                ax_M2_1.plot(sndofs,xI[imotion,:,iorder],color,linewidth=1.0,label=order)
                ax_M2_2.plot(sndofs,yI[imotion,:,iorder],color,linewidth=1.0,label=order)
                ax_M2_3.plot(sndofs,W[ imotion,:,iorder],color,linewidth=1.0,label=order)



ax_M1_1.legend()
ax_M1_2.legend()
ax_M1_3.legend()

ax_M2_1.legend()
ax_M2_2.legend()
ax_M2_3.legend()

#ax_M1_1.set_xlabel("Spatial DOF's")
#ax_M1_2.set_xlabel("Spatial DOF's")
#ax_M1_3.set_xlabel("Spatial DOF's")
#
#ax_M2_1.set_xlabel("Spatial DOF's")
#ax_M2_2.set_xlabel("Spatial DOF's")
#ax_M2_3.set_xlabel("Spatial DOF's")
ax_M1_1.set_xlabel("Mesh Index")
ax_M1_2.set_xlabel("Mesh Index")
ax_M1_3.set_xlabel("Mesh Index")

ax_M2_1.set_xlabel("Mesh Index")
ax_M2_2.set_xlabel("Mesh Index")
ax_M2_3.set_xlabel("Mesh Index")

ax_M1_1.set_ylabel('Impulse-X')
ax_M1_2.set_ylabel('Impulse-Y')
ax_M1_3.set_ylabel('Work')

ax_M2_1.set_ylabel('Impulse-X')
ax_M2_2.set_ylabel('Impulse-Y')
ax_M2_3.set_ylabel('Work')



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

#ax_M1_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M1_4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M1_7.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_5.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_7.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M2_8.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_5.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_7.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax_M3_8.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


#ax_M1_1.set_ylim((-0.1e-06,0.1e-06))
#ax_M1_2.set_ylim((-0.012,-0.011))
#ax_M1_3.set_ylim((-0.009,-0.008))
#ax_M1_4.set_ylim((-0.1e-06,0.1e-06))
#ax_M1_5.set_ylim((-0.014,-0.01))
#ax_M1_6.set_ylim((-0.009,-0.008))
#ax_M1_7.set_ylim((-0.1e-06,0.1e-06))
#ax_M1_8.set_ylim((-0.015,-0.005))
#ax_M1_9.set_ylim((-0.012,-0.005))
#
#ax_M2_1.set_ylim((-0.1e-6,0.1e-6))
#ax_M2_2.set_ylim((-0.1e-6,0.1e-6))
#ax_M2_3.set_ylim((-0.1e-6,0.1e-6))
#ax_M2_4.set_ylim((-0.1e-6,0.1e-6))
#ax_M2_5.set_ylim((-0.1e-6,0.1e-6))
#ax_M2_6.set_ylim((-0.12,-0.1))
#ax_M2_7.set_ylim((-0.1e-6,0.1e-6))
#ax_M2_8.set_ylim((-0.1e-6,0.1e-6))
#ax_M2_9.set_ylim((-0.20,-0.15))
#
#ax_M3_1.set_ylim((-0.1e-06,0.1e-06))
#ax_M3_2.set_ylim((-0.1e-06,0.1e-06))
#ax_M3_3.set_ylim((-1.e-5,-2.e-4))
#ax_M3_4.set_ylim((-0.1e-06,0.1e-06))
#ax_M3_5.set_ylim((-0.1e-06,0.1e-06))
#ax_M3_6.set_ylim((-0.0006,-0.0003))
#ax_M3_7.set_ylim((-0.1e-06,0.1e-06))
#ax_M3_8.set_ylim((-0.1e-06,0.1e-06))
#ax_M3_9.set_ylim((-0.04,-0.03))
#
#ax_M4_1.set_ylim((-0.006,-0.005))
#ax_M4_2.set_ylim((-0.010,-0.008))
#ax_M4_3.set_ylim((-0.006,-0.004))
#ax_M4_4.set_ylim((-0.006,-0.004))
#ax_M4_5.set_ylim((-0.009,-0.006))
#ax_M4_6.set_ylim((-0.130,-0.12))
#ax_M4_7.set_ylim((-0.008,-0.004))
#ax_M4_8.set_ylim((-0.008,-0.004))
#ax_M4_9.set_ylim((-0.20,-0.12))






M1.savefig('UCB_convergence_motion1.png', bbox_inches='tight', dpi=800)
M2.savefig('UCB_convergence_motion2.png', bbox_inches='tight', dpi=800)

plt.show()




