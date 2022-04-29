import numpy as np
import matplotlib.pyplot as plt
import os

# plt.rc("font",family="serif",serif="Computer Modern")
# plt.rc("text",usetex=True)
# plt.rcParams.update({'mathtext.default': 'regular'})

### input parameters
mdm = 0.1
gD = 0.001
rd = 0.12

#mzpl = np.logspace(-3,-0.1,50)
mdml = np.logspace(-3,0,50)
#epsl_sq = np.logspace(-6,-1,50)
epsl_sq = [1e-2]
#mdml = np.logspace(-3,0,50)
#mzp = 0.3
Xfl = []
nChanl = []

#omegal = np.zeros((len(epsl_sq),len(mzpl)))
omegal = np.zeros((len(epsl_sq),len(mdml)))

bad_results = [] #tuples with (mzp,eps,gzp) of bad results (omega=-1)
mzpi = 0
#for mzp in mzpl:
for mdm in mdml:
	mzp = mdm*3
	epsi = 0
	for eps_sq in epsl_sq:
		print("mzp = %.2e, eps2 = %.2e" %(mzp,eps_sq))
		eps = np.sqrt(eps_sq)
		gzp = eps*0.31343
		### write input file
		input_file = "data_python.par"

		Mdm = "Mdm	%.2E" % mdm
		MZp = "MZp	%.2E" % mzp
		gZp = "gZp	%.2E" % gzp
		gZpDm = "gZpDm	%.1E" % gD

		input_parameters = [Mdm,MZp,gZp,gZpDm]
		np.savetxt(input_file,input_parameters, delimiter=" ", fmt="%s")

		### run main.c
		os.system("./mymain " + input_file)

		Xf_file,omega_file = np.loadtxt("Omega.dat")
		if omega_file == -1:
			bad_results.append((mzp,eps_sq,gzp))
		omegal[epsi][mzpi] = omega_file
		Xfl.append(Xf_file)

		### read channels file
		Chan = np.loadtxt("channels/channels_mdm%.2e.dat"%mdm,dtype='str')
		if np.shape(Chan)==(5,):
			nChan = 1
		else:
			nChan = len(Chan)
		nChanl.append(nChan)

		epsi +=1 
	mzpi +=1

### plot contour
# plt.figure()
# ax = plt.gca()
# leveli = np.floor(np.log10(np.min(np.where(omegal > 0, omegal, np.inf))))
# levelf = np.ceil(np.log10(np.max(omegal)))	
# levels=np.power(10,np.arange(leveli,levelf+1,0.1))
# CL_f = plt.contourf(np.log10(mzpl),np.log10(epsl_sq),np.log10(omegal),levels=np.log10(levels))
# CL  = plt.contour(np.log10(mzpl),np.log10(epsl_sq),np.log10(omegal),[np.log10(rd)],colors='k',linestyles='solid')
# cbar=plt.colorbar(CL_f)
# cbarticks = np.arange(leveli,levelf+1,2)
# cbarlabels = ['$10^{%d}$'%(v) for v in np.arange(leveli,levelf+2,2)]
# xticks = np.arange(np.ceil(np.log10(mzpl[0])),np.log10(mzpl[-1])+1,1)
# yticks = np.arange(np.ceil(np.log10(epsl_sq[0])),np.log10(epsl_sq[-1]),1)
# xlabels = ['$10^{%d}$'%(v) for v in xticks]
# ylabels = ['$10^{%d}$'%(v) for v in yticks]
# ax.set_xticks(xticks)
# ax.set_yticks(yticks)
# ax.set_xticklabels(xlabels)
# ax.set_yticklabels(ylabels)
# cbar.set_ticks(cbarticks)
# cbar.set_ticklabels(cbarlabels)
# cbar.ax.set_ylabel(r'$\Omega h^2$',fontsize='12')
# plt.xlabel(r"$m_{Z^\prime}$ (GeV)")
# plt.ylabel(r"$\epsilon^2$")
# plt.text(-2.5,np.log10(min(epsl_sq))+0.5,r"$m_{\chi} = %.1f$ GeV"%mdm)
# plt.text(-2.5,np.log10(min(epsl_sq))+1.0,r"$g_{\chi} = 10^{%d}$"%np.log10(gD))

#cbar.ax.tick_params(labelsize=12)
#plt.savefig("figures/contour_beps1e-6_fast0.pdf")
# for i in yl:
# 	plt.axhline(y=np.log10(i),c='k')
# for i in mdml:
# 	plt.axvline(x=np.log10(i),c='k')

# ### suplot Xf and omega
# plt.figure()
# plt.subplot(2,1,1)
# plt.loglog(mdml,omegal[0],label="Omega")
# plt.legend()
# plt.subplot(2,1,2)
# plt.loglog(mdml,Xfl,label="Xf",c='red')
# plt.legend()


### subplot omega and channels
plt.figure()

plt.subplot(2,1,1)
plt.title(r"$g_D=%.1e$"%gD)
plt.loglog(mdml,omegal[0])
plt.ylabel(r'$\Omega h^2$')
plt.subplot(2,1,2)
plt.semilogx(mdml,nChanl,label="N of channels",c='k')
plt.ylabel("number of channels")
plt.xlabel(r'$m_\chi$ (GeV)')



plt.show()


