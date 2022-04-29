import numpy as np
import matplotlib.pyplot as plt
import os
import sys

######## Run Madgraph
run_name = "DarkPhoton_Ebeam4_SF1"
os.chdir("/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0")
MG_dir = "./bin/mg5_aMC"
input_dir = "/Users/Taylor/Documents/PhD/Code/MGscripts/input_file.txt"

### run settings
scale_factor = 1. # multiplies all relevant energy scales
N_event = 5000
ebeam1 = 4.0*scale_factor
ebeam2 = 171.0 *scale_factor

mAp_list = [0.001,0.01,0.1,1.0,2.0]
#mAp_list = [10**(i/10) for i in range(-30,10,1)]
#mAp_list = [0.01]
mAp_list = [x * scale_factor for x in mAp_list]

Mnt = 171.0 * scale_factor
Me = 0.000511 * scale_factor

### write input file
l1 = "import model DarkPhoton_UFO"
l2 = "generate e- nt > e- nt Ap / z h"
l3 = "output "+run_name
l4 = "launch"
l5 = "set nevents = %d"%N_event
l6 = "set ebeam1 %.3f"%ebeam1
l7 = "set ebeam2 %.3f"%ebeam2
l8 ="set use_syst False"
l9 = "set lpp1 0"
l10 = "set lpp2 0"
l11 = "set ptl -1"
l12 = "set etal -1"
l13 = "set MAp scan:"+str(mAp_list)
l14 = "set WAp auto"
l15 = "set Mnt %.3f"%Mnt
l16 = "set Me %.3f"%Me
l17 = "set ptl 0.01"
# l16 = "set {2001:2.5} = eta_min_pdg"
# l17 = "set {2001:2.5} = eta_max_pdg"

input_parameters = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16]
np.savetxt(input_dir,input_parameters, delimiter=" ", fmt="%s")

### run MG
os.system(MG_dir +" "+ input_dir)

os.chdir("/Users/Taylor/Documents/PhD/Code/MGscripts")
