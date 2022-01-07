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


M1 = plt.figure(figsize=(5,4))
M1.set_facecolor('White')
ax = M1.add_subplot(111)


# AFRL Cases 
afrl_grids   = ['h2']
afrl_motions = ['M1']
afrl_physics = ['Re10','Euler']
afrl_orders  = ['p1', 'p2', 'p3']
afrl_times   = ['0.001']

for imotion,motion in zip(range(len(afrl_motions)),afrl_motions):
    for iphysic,physic in zip(range(len(afrl_physics)),afrl_physics):


        if physic == 'Re10':
            motion1_truth = ['h2','p3','0.001']
            motion2_truth = ['h2','p3','0.001']
            motion3_truth = ['h2','p3','0.001']
            motion4_truth = ['h2','p3','0.001']
        elif physic == 'Re1000':
            motion1_truth = ['h2','p3','0.001']
            motion2_truth = ['h2','p3','0.001']
            motion3_truth = ['h2','p3','0.001']
            motion4_truth = ['h2','p3','0.001']
        elif physic == 'Euler':
            motion1_truth = ['h1','p3','0.001']
            motion2_truth = ['h1','p3','0.001']
            motion3_truth = ['h1','p3','0.001']
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

        if (physic == 'Euler'):
            afrl_file = 'AFRL/'+motion+'-'+"ReInf"+'-'+h_ref+'-'+p_ref+'/'+'cyl-t'+t_ref
        else:
            afrl_file = 'AFRL/'+motion+'-'+physic+'-'+h_ref+'-'+p_ref+'/'+'cyl-t'+t_ref
    
        # Load data
        if (os.path.isfile(afrl_file)):
            afrl_data = np.loadtxt(afrl_file, skiprows=1)
        else:
            print('Data not found: '+afrl_file)
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

            if physic == "Euler":
                color = 'b--'
            elif physic == "Re1000":
                color = 'r--'
            elif physic == "Re10":
                color = 'g--'
            
            if (motion == 'M1'):
                ax.plot(afrl_x,afrl_Fy,  color,linewidth=1.0,label=physic)


t = np.linspace(0,2.,100)
density = 1.
volume  = np.pi * 0.5 * 0.5
mass = density*volume
accel = 3.*t - (9./4.)*t*t
F = - mass * accel
ax.plot(t,F, 'k', linewidth=1.0, label="F=ma")



ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Force-Y')



#M1.set_xlim((0.,2.))
#M1.set_ylim((-0.1e-03,0.1e-03))


M1.savefig('ForceComparison.png', bbox_inches='tight', dpi=800)

plt.show()




