### Plot Histogram of Momentums from madgraph output LHE files

import sys
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.integrate as spyint

run_name = "10GeV_DP"

def event_file_dir(run_name,run_str):
	return "/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0/"+run_name+"/Events/run_"+run_str+"/unweighted_events.lhe.gz"

nbins = 50

particles = [r"$e^-_{i}$",r"$N_i$",r"$e^-_{f}$",r"$N_f$",r"$A^\prime$"] # particles in the process
particle = 2 # index of particle you want to analyze
dist = 'E' # E, pT, px, py, pz

mApl_plot = [0.01,0.1,0.2,0.5,1.0,1.5] # list of mAp to plot
mApl_plot = [1.5] # list of mAp to plot

mApl = mApl_plot #list of mAp that were run
mApl = [1.5,5.0]


colors = ["tab:blue","yellow","tab:green","tab:red","tab:purple","tab:brown"]

sys.path.append('../../')
os.chdir("... MG5_aMC_v3_2_0/madgraph/various") #insert directory to the various folder
from lhe_parser import EventFile, FourMomentum

j=0
plt.figure()
for mAp in mApl_plot:
	if min(mApl) <= mAp <= max(mApl):
		run = np.argmin(abs(np.array(np.log10(mApl)) - np.log10(mAp)))+1
		mAp = mApl[run-1]
	else:
		print(r"$m_A^\prime = %.2e GeV$ was not run." %mAp)

	if run < 9 :
		run_str = "0%d"%run
	else:
		run_str = "%d"%run

	lhe = EventFile(event_file_dir(run_name,run_str))
	i=0
	data = []
	for event in lhe:
		p = [FourMomentum(particle) for particle in event]
		if dist == 'E':
			value = p[particle].E # E - energy
			plt.xlabel(r"$E$ (GeV)")
		if dist == 'pT':
			value = np.sqrt(abs(p[particle].px**2+p[particle].py**2)) # pT - transverse momentum
			plt.xlabel(r"$|p_T|$ (MeV)")
		if dist == 'pz':
			value = p[particle].pz
			plt.xlabel(r"$p_z$ (GeV)")
		if dist == 'px':
			value = abs(p[particle].px)
			plt.xlabel(r"$p_x$ (GeV)")
		data.append(value)
		i+=1
	data = np.array(data)
	range_i = (min(data),max(data))
	mAp_MeV = mAp*1000

	if dist == "pT":
		zord = 10-j
	else:
		zord = j

	#plot histogram
	plt.hist(data,bins=nbins,color=colors[j],histtype="step",density=True,range=range_i,zorder=zord) #outline
	plt.hist(data,bins=nbins,label=r"$m_{A^\prime} = %d$ MeV"%mAp_MeV,color=colors[j],histtype="stepfilled",density=True,alpha=0.3,range=range_i,zorder=zord)
	j+=1

plt.ylabel("counts")
plt.legend()
ax = plt.gca()
ax.set_yscale('log')
plt.title(particles[particle]+" distribution")
os.chdir("/Users/Taylor/Documents/PhD/Code/MGscripts")
#plt.savefig("histograms/E_efinal_mApscan.pdf")
plt.show()
