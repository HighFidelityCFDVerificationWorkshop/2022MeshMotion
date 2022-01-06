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


M1 = plt.figure(figsize=(10,3))
M1.set_facecolor('White')
ax_M1_1 = M1.add_subplot(131)
ax_M1_2 = M1.add_subplot(132)
ax_M1_3 = M1.add_subplot(133)

M2 = plt.figure(figsize=(10,3))
M2.set_facecolor('White')
ax_M2_1 = M2.add_subplot(131)
ax_M2_2 = M2.add_subplot(132)
ax_M2_3 = M2.add_subplot(133)



# UM Cases 
case_grids          = ['h0', 'h1', 'h2', 'h3']
case_grid_sizes     = [520, 1963, 7781, 15882]
case_motions        = ['M1', 'M2']
case_orders         = ['p1', 'p2', 'p3']
case_order_integers = [1, 2, 3]
case_times          = ['naca0', 'naca1', 'naca2', 'naca3']

# Storage (motions,physics,h,p,t)
xI = np.zeros((len(case_motions),len(case_grids),len(case_orders),len(case_times)))
yI = np.zeros((len(case_motions),len(case_grids),len(case_orders),len(case_times)))
zI = np.zeros((len(case_motions),len(case_grids),len(case_orders),len(case_times)))
W  = np.zeros((len(case_motions),len(case_grids),len(case_orders),len(case_times)))
m  = np.zeros((len(case_motions),len(case_grids),len(case_orders),len(case_times)))

for imotion,motion in zip(range(len(case_motions)),case_motions):
    for igrid,grid in zip(range(len(case_grids)),case_grids):
        for iorder,order in zip(range(len(case_orders)),case_orders):
            for itime,time in zip(range(len(case_times)),case_times):

                case_file = 'Tri-'+motion+'-'+grid+'-'+order+'/'+time+'_TimeHist.txt'

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

                    xI[imotion,igrid,iorder,itime] = integrate.simps(case_Fx,  x=case_t)
                    yI[imotion,igrid,iorder,itime] = integrate.simps(case_Fy,  x=case_t)
                    W[ imotion,igrid,iorder,itime] = integrate.simps(case_Wint,x=case_t)




itime = 3
for imotion,motion in zip(range(len(case_motions)),case_motions):
    for iorder,order in zip(range(len(case_orders)),case_orders):

        if order == 'p1':
            sndofs = [1560,5889,23343,47646]
        elif grid == 'p2':
            sndofs = [3120,11778,46686,95292]
        elif grid == 'p3':
            sndofs = [5200,19630,77810,158820]

            
        if order == 'p1':
            color = 'ro-'
        elif order == 'p2':
            color = 'bo-'
        elif order == 'p3':
            color = 'go-'
        elif order == 'p4':
            color = 'ko-'
    
        if (motion == 'M1'):
                ax_M1_1.semilogx(sndofs,xI[imotion,:,iorder,itime],color,linewidth=1.0,label=order)
                ax_M1_2.semilogx(sndofs,yI[imotion,:,iorder,itime],color,linewidth=1.0,label=order)
                ax_M1_3.semilogx(sndofs,W[ imotion,:,iorder,itime],color,linewidth=1.0,label=order)
        elif (motion == 'M2'):
                ax_M2_1.semilogx(sndofs,xI[imotion,:,iorder,itime],color,linewidth=1.0,label=order)
                ax_M2_2.semilogx(sndofs,yI[imotion,:,iorder,itime],color,linewidth=1.0,label=order)
                ax_M2_3.semilogx(sndofs,W[ imotion,:,iorder,itime],color,linewidth=1.0,label=order)



ax_M1_1.legend()
ax_M1_2.legend()
ax_M1_3.legend()

ax_M2_1.legend()
ax_M2_2.legend()
ax_M2_3.legend()

ax_M1_1.set_xlabel("Spatial DOF's")
ax_M1_2.set_xlabel("Spatial DOF's")
ax_M1_3.set_xlabel("Spatial DOF's")

ax_M2_1.set_xlabel("Spatial DOF's")
ax_M2_2.set_xlabel("Spatial DOF's")
ax_M2_3.set_xlabel("Spatial DOF's")


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


#ax_M1_1.set_ylim((-0.1e-08,0.1e-08))
#ax_M1_2.set_ylim((-0.012,-0.011))
#ax_M1_3.set_ylim((-0.009,-0.008))
#
#ax_M2_1.set_ylim((-0.1e-8,0.1e-8))
#ax_M2_2.set_ylim((-0.1e-8,0.1e-8))
#ax_M2_3.set_ylim((-0.1e-8,0.1e-8))


M1.savefig('UM_convergence_motion1.png', bbox_inches='tight', dpi=800)
M2.savefig('UM_convergence_motion2.png', bbox_inches='tight', dpi=800)

plt.show()




