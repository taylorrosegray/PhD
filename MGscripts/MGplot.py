import sys
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.integrate as spyint

run_name = "DarkPhoton_Ebeam4_SF1"
cross_section = 1
kinematic = 0

me = 0.000511 #GeV electron mass
E0 = 4.0 #GeV electron beam energy 
alpha = 1/137
eps = 1.0
scale_factor = 1.
nbins = 50

### list of mAp to plot
mAp_list = [0.001,0.01,0.1,1.0,2.0]
#mAp_list = [10**(i/10) for i in range(-30,10,1)]
mAp_list = [x * scale_factor for x in mAp_list]

if kinematic == 1:
	### Plot Histogram of Momentums

	particle = 2
	particles = [r"$e^-_{i}$",r"$N_i$",r"$e^-_{f}$",r"$N_f$",r"$A^\prime$"]
	dist = 'E' # E, pT, px, py, pz

	mApl = mApl_plot

	#mApl = [10**(i/10) for i in range(-30,10,1)] # list of Ap masses that were run
	#mApl = [0.01,0.1,0.2,0.5,1.0,1.5]
	#mApl = [x * scale_factor for x in mApl]
	#mApl = mApl_plot

	if len(mApl) > 1:
		colors = ["tab:blue","yellow","tab:green","tab:red","tab:purple","tab:brown"]
	else:
		colors = ["purple"]

	sys.path.append('../../')
	os.chdir("/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0/madgraph/various")
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

		lhe = EventFile("/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0/"+run_name+"/Events/run_"+run_str+"/unweighted_events.lhe.gz")
		i=0
		#data = np.zeros((len(lhe)))
		data = []
		for event in lhe:
			p = [FourMomentum(particle) for particle in event]
			if dist == 'E':
				value = p[particle].E # E - energy
				plt.xlabel(r"$E$ (GeV)")
			if dist == 'pT':
				#if p[particle].E < 1.2:
				value = np.sqrt(abs(p[particle].px**2+p[particle].py**2)) # pT - transverse momentum
				plt.xlabel(r"$|p_T|$ (MeV)")
				# else:
				# 	continue
			if dist == 'pz':
				value = p[particle].pz
				plt.xlabel(r"$p_z$ (GeV)")
			if dist == 'px':
				value = abs(p[particle].px)
				plt.xlabel(r"$p_x$ (GeV)")

			#data[i] = value
			data.append(value)
			i+=1
		data = np.array(data)
		data = data/scale_factor
		if dist=="E":
			#range_i = (min(data),max(4,max(data)))
			range_i = (min(data),max(data))
		elif dist == "pT":
			data = data*1000 # convert to MeV
			#range_i = (min(data),800)
			range_i = (min(data),max(data))
		else:
			#range_i = (min(data)-0.1,max(data)+0.1)
			range_i = (min(data),max(data))
		mAp_MeV = mAp*1000/scale_factor

		if dist == "pT":
			zord = 10-j
		else:
			zord = j
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
	# plt.ylim(8e-4,1)
	plt.show()

if cross_section == 1:
	## Plot mass scan
	GV = 1.0 #dark coupling
	EOT = 4e14 #electrons on target
	X0 = 6.76 #radiation lengths g/cm^2 of tungsten
	T = 0.1
	A = 183.84 #atomic mass g/mol
	Ne = 4e14 #electrons on target
	N0 = 6.02e23 #avogadro
	conversion = 1e-6 # 1e-6 micro = 1 p

	run_length = "%d"%len(mAp_list)
	if len(mAp_list) < 9 :
		run_str = "0[1-"+run_length+"]"
	else:
		run_str = "[01-"+run_length+"]"

	path = "/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0/"+run_name+"/Events/scan_run_"+run_str+".txt"
	mAp, cross_section, width = np.loadtxt(path,skiprows=1,usecols=(1,2,3),unpack = True) # GeV, pb
	mAp = mAp/scale_factor
	events = cross_section*1e-36* T * X0 * Ne * N0 / A
	cross_section = cross_section*conversion*scale_factor**2
	#mAp_paper, cross_section_paper = np.loadtxt("cs_vs_mAp.txt",unpack=True) # MeV, events/eps^2

	plt.figure()
	plt.rcParams["mathtext.fontset"] = "cm"
	#plt.loglog(mAp,X0*EOT*N0*conversion/A*cross_section*0.000001 / GV**2,label="MadGraph") #convert cross section to micro b
	#plt.loglog(mAp,cross_section,label="MG")
	#plt.loglog(mAp_paper/1000,cross_section_paper,label="LDMX paper")
	plt.loglog(mAp*1000,events)
	plt.xlabel(r"$m_{A^{\prime}}$ (MeV)")
	#plt.ylabel(r"$\sigma (\mu b)$")
	plt.ylabel("Events")
	

	# ### theoretical cross section
	# def G2(t): #form factor
	# 	Z = 74
	# 	A = 184
	# 	a = 111 * Z**(1/3)/me
	# 	d = 0.164 * A**(-2/3)
	# 	return (a**2 * t / (1+a**2*t))**2 * (1/(1+t/d))**2 * Z**2

	# def cross_section_th(mAp): # theoretical A' cross section
	# 	a = max(me**2/mAp**2,mAp**2/E0**2)
	# 	beta = np.sqrt(1-mAp**2/E0**2)
	# 	tmin = 0.0
	# 	def chi_int(t):
	# 		return G2(t)*(t-tmin)/t**2
	# 	chi = spyint.quad(chi_int,0, mAp**2)[0]
	# 	result = 8/3 * alpha**3 * eps**2 * beta / mAp**2 * chi * np.log(1/a)
	# 	return result
	# csl = []
	# for m in mAp:
	# 	cs = cross_section_th(m)*3.89379e8
	# 	csl.append(cs*1e-6) #convert from pb to micro-b
	# plt.loglog(mAp,csl,label="theoretical")
	# #print(r"sigma_th(mAp = %.4f) = "%mAp,cs, r'$pb$')

	plt.legend()
	plt.show()
