import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

### input parameters
alphad = 0.5
gD_value = np.sqrt(4*np.pi*alphad)
N_points = 50
m_ratio = 3.0 # mAp/mdm
b5 = 0.0
b6R = 0.0
b6I = 0.0
b7R = 0.0
b7I = gD_value
run_name = "vectorDM_b7I_ratio"+str(m_ratio)+"_alphaD"+str(alphad)

### loop over
mdml = np.logspace(-3,0.5,N_points)
yl = np.logspace(-16,-7,N_points)

### initialize
omegal = np.zeros((len(yl),len(mdml)))
Xfl = []
nChanl = []

rd = 0.12
ee = 0.31343
input_file = "data_python.par"

def gSM_f(y,alphad,m_ratio):
	return np.sqrt(y*ee**2*m_ratio**4/alphad)

###### loop and run micromegas
mdmi = 0
for mdm in mdml:
	mAp = mdm*m_ratio
	yi = 0
	for y in yl:
		print("mdm = %.2e, y = %.2e" %(mdm,y))
		gSM = gSM_f(y,alphad,m_ratio)

		### write input file

		l1   = "mdm    %.2E" % mdm
		l2   = "mAp    %.2E" % mAp
		l3   = "gSM    %.3E" % gSM
		l4   = "b5     %.2E" %b5
		l5   = "b6R    %.2E" %b6R
		l6   = "b6I    %.2E" %b6I
		l7   = "b7R    %.2E" %b7R
		l8   = "b7I    %.2E" %b7I


		input_parameters = [l1,l2,l3,l4,l5,l6,l7,l8]
		np.savetxt(input_file,input_parameters, delimiter=" ", fmt="%s")

		### run main.c
		os.system("./mymain " + input_file)

		### read omega file
		Xf,omega_file = np.loadtxt("Omega.dat",unpack=True)
		#omegal.append(omega_file)
		omegal[yi][mdmi] = omega_file
		Xfl.append(Xf)

		### read sigmav file
		# T_file,sigmav_file = np.loadtxt("sigmav.dat",unpack=True)
		# plt.figure()
		# plt.loglog(T_file,sigmav_file)
		# plt.title(r"$m_\chi = %.2e$, y = %.2e" % (mdm,y))
		# plt.savefig("figures/sigmav/mdm%.2e_y%.2e.jpeg"% (mdm,y))
		# plt.close()

		### read channels file
		# Chan = np.loadtxt("channels/channels_mdm%.2e_gSM%.3e.dat"%(mdm,gSM),dtype='str')
		# if np.shape(Chan)==(5,):
		# 	nChan = 1
		# else:
		# 	nChan = len(Chan)
		# nChanl.append(nChan)

		yi += 1
	mdmi +=1

### save data
path = "runs/"+run_name+"/"
if not os.path.isdir(path):
	os.mkdir(path)
np.save(path+"omega",omegal)
np.save(path+"y",yl)
np.save(path+"mdm",mdml)

#### Plotting
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.rc("font",family="serif",serif="Computer Modern")
plt.rc("text",usetex=True)
plt.rcParams.update({'mathtext.default': 'regular'})

if np.all(np.isnan(omegal)) or np.all(omegal == -1.0):
	print("ALL NANS :(")
else:
	### plot contour
	plt.figure()
	ax = plt.gca()
	leveli = np.floor(np.log10(np.min(np.where(omegal > 0, omegal, np.inf))))
	levelf = np.ceil(np.log10(np.max(omegal)))	
	levels=np.power(10,np.arange(leveli,levelf+1,0.1))
	CL_f = plt.contourf(np.log10(mdml),np.log10(yl),np.log10(omegal),levels=np.log10(levels))
	CL  = plt.contour(np.log10(mdml),np.log10(yl),np.log10(omegal),[np.log10(rd)],colors='k',linestyles='solid')
	cbar=plt.colorbar(CL_f)
	cbarticks = np.arange(leveli,levelf+1,2)
	cbarlabels = ['$10^{%d}$'%(v) for v in np.arange(leveli,levelf+2,2)]
	xticks = np.arange(np.ceil(np.log10(mdml[0])),np.log10(mdml[-1]),1)
	yticks = np.arange(np.ceil(np.log10(yl[0])),np.log10(yl[-1]),1)
	xlabels = ['$10^{%d}$'%(v) for v in xticks]
	ylabels = ['$10^{%d}$'%(v) for v in yticks]
	ax.set_xticks(xticks)
	ax.set_yticks(yticks)
	ax.set_xticklabels(xlabels)
	ax.set_yticklabels(ylabels)
	cbar.set_ticks(cbarticks)
	cbar.set_ticklabels(cbarlabels)
	cbar.ax.set_ylabel(r'$\Omega h^2$',fontsize='12')
	plt.xlabel(r"$m_\chi$ (GeV)")
	plt.ylabel(r"$y$")


	# #add elastic scalar, pseudo-dirac
	mdm_scalar,y_scalar = np.loadtxt("digitized_datasets/scalarelastic.txt",unpack=True,delimiter=",")
	mdm_dirac,y_dirac = np.loadtxt("digitized_datasets/pseudodirac.txt",unpack=True,delimiter=",")
	mdm_LDMX,y_LDMX = np.loadtxt("digitized_datasets/LDMX_projection.txt",unpack=True,delimiter=",")
	# plt.plot(np.log10(mdm_scalar/1000),np.log10(y_scalar),label=r"Scalar Elastic DM",color='black',ls="dashed")
	# plt.plot(np.log10(mdm_dirac/1000),np.log10(y_dirac),label=r"Pseudo-Dirac DM",color='black',ls="dotted")
	# plt.plot(np.log10(mdm_LDMX/1000),np.log10(y_LDMX),label="LDMX Projection",c='tab:pink')

	#legend
	contour_legend_lines = [Line2D([0], [0], color='black'),Line2D([0], [0], color='black',ls='dashed'),Line2D([0], [0], color='black',ls='dotted'),Line2D([0], [0], color='tab:pink')]
	contour_legend_labels = ["Vector DM","Scalar Elastic DM", "Pseudo-Dirac DM", "LDMX Projection"]
	plt.legend(contour_legend_lines,contour_legend_labels,loc="lower right")

	#cbar.ax.tick_params(labelsize=12)
	#plt.savefig("figures/contour_beps1e-6_fast0.pdf")
	# for i in yl:
	# 	plt.axhline(y=np.log10(i),c='k')
	# for i in mdml:
	# 	plt.axvline(x=np.log10(i),c='k')

	plt.text(-0.3,-13.5,r"$\alpha_D = %.1f$"%alphad)
	plt.ylim(-16,-7)
	plt.xlim(-3,np.log10(max(mdm_dirac/1000)))
	plt.title(r"$m_{A^\prime}/m_\chi=%.1f$"%m_ratio)
	plt.savefig(path+"VectorDM.pdf")

	# plt.axvline(x=np.log10(0.00255),c='k')
	# plt.axvline(x=np.log10(2*0.00255),c='k',ls='dashed')

	plt.show()

# plt.figure()
# plt.loglog(mdml,omegal[0])


# plt.subplot(2,1,1)
# plt.title("y=%.1e"%y)
# plt.loglog(mdml,omegal[0])
# plt.ylabel(r'$\Omega h^2$')
# #plt.legend()
# plt.subplot(2,1,2)
# #plt.loglog(mdml,np.array(mdml)/np.array(Xfl),label="Xf",c='red')
# plt.semilogx(mdml,nChanl,label="# channels",c='k')
# plt.ylabel("# channels")
# plt.xlabel(r'$m_\chi$ (GeV)')
# # plt.axvline(x=0.214,label=r"Drop in $g_{eff}$",c='m')
# # plt.axvline(x=0.1,c='m')
# # plt.fill_betweenx([min(nChanl),max(nChanl)],0.1,0.214,color='m',alpha=0.5)
# #plt.legend()


# plt.figure()
# plt.subplot(2,1,1)
# plt.title(r"$m_{\chi}=%.1e$"%mdm)
# plt.loglog(yl, omegal[:,0])
# plt.ylabel(r'$\Omega h^2$')
# plt.subplot(2,1,2)
# plt.semilogx(yl,nChanl,label="# channels",c='k')
# plt.ylabel("# channels")
# plt.xlabel(r'$y$')

### plot gEff
# plt.figure()
# T_file,geff_file = np.loadtxt("gEff.dat",unpack=True)
# plt.loglog(T_file, geff_file)



#plt.show()

# mfe = 0.000511 #mass of electron
# mfm = 0.1066 #mass of muon
# mft = 1.777 #mass of tau
# mfd = 0.007 #mass of down quark
# mfu = 0.003 #mass of up quark
# mfs = 0.120 #mass of strange quark
# mfc = 1.2 #mass of charm quark
# mfb = 4.3 #mass of bottom quark
# mfto = 174 #mass of top quark
# mfl = [mfe,mfm,mft,mfd,mfu,mfs,mfc,mfb,mfto]

