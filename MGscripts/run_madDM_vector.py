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
gD = np.sqrt(4*np.pi*0.5)
b5 = gD
b6I = 0.0
b6R = 0.0
b7I = 0.0
b7R = 0.0
run_name = "VectorDM_madDM_b5"

N_points = 30
ymin = 1e-16
ymax = 1e-5
min_mdm = 1e-3 # GeV
max_mdm = 10.0  # GeV
ratio_mAp_mdm = 3.0

def gSM_f(y):
	return np.sqrt(ee**2 * 4 * np.pi * (ratio_mAp_mdm)**4 / gD**2)*np.sqrt(y)

min_gSM = gSM_f(ymin)
max_gSM = gSM_f(ymax)

# min_gSM = 1e-8
# max_gSM = 1.0
gSMl = np.array([10**(np.log10(min_gSM)+i*(np.log10(max_gSM)-np.log10(min_gSM))/N_points) for i in range(0,N_points,1)])
yl = np.array([10**(np.log10(ymin)+i*(np.log10(ymax)-np.log10(ymin))/N_points) for i in range(0,N_points,1)])


min_mAp = min_mdm*ratio_mAp_mdm
max_mAp = max_mdm*ratio_mAp_mdm


mdml = np.logspace(np.log10(min_mdm),np.log10(max_mdm),N_points)
Omega = np.zeros((N_points,N_points))


### reproduce fig2
# mdml=[20]
# min_mAp = 20.0
# max_mAp = 140.0
# N_points = 200
# b5 = 0.0
# b6I = gD
# b6R = 0.0
# b7I = 0.0
# b7R = 0.0
for j in range(len(mdml)):

	mdm = mdml[j]
	mAp = mdm*ratio_mAp_mdm

	### write input file
	l1 = "import model VectorDM_UFO"
	l2 = "define darkmatter xv"
	l3 = "generate relic_density"
	l4 = "output "+run_name+"/mdm%.3e"%mdm
	l5 = "launch "+run_name+"/mdm%.3e"%mdm
	l6 = "set mdm %.3e"%mdm
	l7 = "set mAp %.3e"%mAp
	#l6 = "set mdm scan:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_mdm),np.log10(max_mdm),np.log10(min_mdm),N_points,N_points) 
	#l7 = "set mAp scan1:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_mAp),np.log10(max_mAp),np.log10(min_mAp),N_points,N_points) 
	l8 = "set gSM scan:[10**(%.3e+i*(%.3e-%.3e)/%d) for i in range(0,%d,1)]"%(np.log10(min_gSM),np.log10(max_gSM),np.log10(min_gSM),N_points,N_points) 
	#l8 = "set gSM 0.5"
	l9 = "set b5 %.3e"%b5
	l10 = "set b6R %.3e"%b6R
	l11 = "set b6im %.3e"%b6I
	l12 = "set b7R %.3e"%b7R
	l13 = "set b7im %.3e"%b7I

	input_parameters = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13]
	np.savetxt(input_dir,input_parameters, delimiter=" ", fmt="%s")


### run madDM

	os.system("python2 "+madDM_dir +" "+ input_dir)

os.chdir("/Users/Taylor/Documents/PhD/Code/MGscripts")



###taacas_check(min_x, max_x, nsteps, logscale):