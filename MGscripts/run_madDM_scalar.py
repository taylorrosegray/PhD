import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.lines import Line2D

## variables
ee = 0.31343

##### run madDM to compute relic abundance for VectorDM_UFO

os.chdir("/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0")
madDM_dir = "bin/maddm.py"
input_dir = "/Users/Taylor/Documents/PhD/Code/MGscripts/input_file.txt"

##### run settings

alphaD = 0.5
gD = np.sqrt(4*np.pi*alphaD)
#gSM = 0.5
run_name = "ScalarDM_madDM_DMsimp_s_spin1_testing"

N_points = 15
ymin = 1e-16
ymax = 1e-2
min_mdm = 1e-3 # GeV
max_mdm = 10.0  # GeV
ratio_mAp_mdm = 3.0
#ratio_mAp_mdm = 10.0


### defines masses/couplings to loop over
def gSM_f(y):
	return np.sqrt(ee**2 * 4 * np.pi * (ratio_mAp_mdm)**4 / gD**2)*np.sqrt(y)

min_gSM = gSM_f(ymin)
max_gSM = gSM_f(ymax)
gSMl = np.array([10**(np.log10(min_gSM)+i*(np.log10(max_gSM)-np.log10(min_gSM))/N_points) for i in range(0,N_points,1)])
yl = np.array([10**(np.log10(ymin)+i*(np.log10(ymax)-np.log10(ymin))/N_points) for i in range(0,N_points,1)])

min_mAp = min_mdm*ratio_mAp_mdm
max_mAp = max_mdm*ratio_mAp_mdm

mdml = np.logspace(np.log10(min_mdm),np.log10(max_mdm),N_points)
Omega = np.zeros((N_points,N_points))

mdml = [0.001,0.01,0.1,1.0]
for j in range(len(mdml)):

	mdm = mdml[j]
	mAp = mdm*ratio_mAp_mdm

	### write input file
	# l1 = "import model ScalarDM_UFO"
	# l2 = "define darkmatter xc"
	# l3 = "generate relic_density"
	# l4 = "output "+run_name+"/mdm%.3e"%mdm
	# l5 = "launch "+run_name+"/mdm%.3e"%mdm
	# l6 = "set mdm %.3e"%mdm
	# l7 = "set mAp %.3e"%mAp
	# #l6 = "set mdm scan:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_mdm),np.log10(max_mdm),np.log10(min_mdm),N_points,N_points) 
	# #l7 = "set mAp scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_mAp),np.log10(max_mAp),np.log10(min_mAp),N_points,N_points) 
	# #l8 = "set gSM scan:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points) 
	# l8 = "set gSM 0.01"
	# l9 = "set gD %.3e"%gD

	### DMsimp_s_spin1
	l1 = "import model DMsimp_s_spin1"
	l2 = "define darkmatter xc"
	l3 = "generate relic_density"
	l4 = "output "+run_name+"/mdm%.3e"%mdm
	l5 = "launch "+run_name+"/mdm%.3e"%mdm
	l6 = "set MXc %.3e"%mdm
	l7 = "set MY1 %.3e"%mAp
	#l8 = "set gSM scan:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points) 
	l8 = "set gAd11 0.0"
	l9 = "set gAu11 0.0"
	l10 = "set gAd22 0.0"
	l11 = "set gAu22 0.0"
	l12 = "set gAd33 0.0"
	l13 = "set gAu33 0.0"
	l14 = "set gAu31 0.0"
	l15 = "set gAd31 0.0"

	l22 = "set gVl11 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	#l22 = "set gVl11 1.0"

	l16 = "set gVd11 0.0"
	l17 = "set gVu11 0.0"
	l18 = "set gVd22 0.0"
	l19 = "set gVu22 0.0"
	l20 = "set gVd33 0.0"
	l21 = "set gVu33 0.0"
	l23 = "set gVl22 0.0"
	l24 = "set gVl33 0.0"


	# l16 = "set gVd11 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	# l17 = "set gVu11 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	# l18 = "set gVd22 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	# l19 = "set gVu22 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	# l20 = "set gVd33 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	# l21 = "set gVu33 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	# l23 = "set gVl22 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)
	# l24 = "set gVl33 scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points)


	l25 = "set gVu31 0.0"
	l26 = "set gVd31 0.0"
	l27 = "set gAl11 0.0"
	l28 = "set gAl22 0.0"
	l29 = "set gAl33 0.0"
	l30 = "set gnu11 0.0"
	l31 = "set gnu22 0.0"
	l32 = "set gnu33 0.0"
	l33 = "set gVXc %.3e"%(gD*2)
	#l33 = "set gVXc scan:[1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1.0]"
	l34 = "set gVXd 0.0"
	l35 = "set gAXd 0.0"
	l36 = "set gVu31 0.0"
	l37 = "set gVd31 0.0"

	input_parameters = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37]
	np.savetxt(input_dir,input_parameters, delimiter=" ", fmt="%s")

### run madDM

	os.system("python2 "+madDM_dir +" "+ input_dir)

os.chdir("/Users/Taylor/Documents/PhD/Code/MGscripts")

###taacas_check(min_x, max_x, nsteps, logscale):