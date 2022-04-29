import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.lines import Line2D
ee = 0.31343

###### analyze results ######

##### run settings

run_name = "ScalarDM_madDM_DMsimp_s_spin1_testing"

N_points = 4
alphaD = 0.5

gD = np.sqrt(4*np.pi*0.5)
ymin = 1e-16
#ymax = 1e-5
ymax = 1e-2
min_mdm = 1e-3 # GeV
max_mdm = 10.0  # GeV
ratio_mAp_mdm = 3.0
mdml = np.logspace(np.log10(min_mdm),np.log10(max_mdm),N_points)
mdml = [0.001,0.01,0.1,1.0]

def gSM_f(y):
	return np.sqrt(ee**2 * 4 * np.pi * (ratio_mAp_mdm)**4 / gD**2)*np.sqrt(y)

min_gSM = gSM_f(ymin)
max_gSM = gSM_f(ymax)
gSMl = np.array([10**(np.log10(min_gSM)+i*(np.log10(max_gSM)-np.log10(min_gSM))/N_points) for i in range(0,N_points,1)])
yl = np.array([10**(np.log10(ymin)+i*(np.log10(ymax)-np.log10(ymin))/N_points) for i in range(0,N_points,1)])

Omega = np.zeros((N_points,N_points))
sigmav_xf = np.zeros((N_points,N_points))
xf = np.zeros((N_points,N_points))

###single parameter
#results_path = "/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0/"+run_name+"/output/run_01/MadDM_results.txt" 
#results = np.loadtxt(results_path,delimiter="=",skiprows=9,usecols=1,dtype=str)
#RD = float(results[0]) # relic density
#print("relic density = %.2e"%RD)
plt.figure()
colors=['r','b','g','orange']

##parameter scan
for j in range(N_points):
	results_path = "/Users/Taylor/Documents/PhD/Code/MG5_aMC_v3_2_0/"+run_name+"/mdm%.3e/output/scan_run_01.txt"% mdml[j] 
	results_param,results_Om,results_xf,results_sigmavxf = np.loadtxt(results_path,skiprows=15,usecols=[1,2,3,4],unpack=True)
	#results_param,results_Om,results_xf,results_sigmavxf = np.loadtxt(results_path,skiprows=15,usecols=[1,10,11,12],unpack=True)
	#results_mu,results_ap,results_s,results_d,results_c,results_tau,results_e,results_t,results_b,results_u = np.loadtxt(results_path,skiprows=15,usecols=[5,6,7,8,9,10,11,12,13,14],unpack=True) #contribution of channels
	# if N_points > 1:
	# 	Omega[:,j] = results_Om # access Omega by Omega[y_index][mdm_index]
	# 	sigmav_xf[:,j] = results_sigmavxf
	# 	xf[:,j] = results_xf
	#### plot relic density vs parameter
	ax = plt.gca()
	plt.loglog(results_param,results_Om,c=colors[j],label="mdm=%.3e"%mdml[j])
	
	#plt.xlabel(r"$m_{A^\prime}$ (GeV)")
	plt.xlabel("gSM")
	plt.ylabel(r"$\Omega h^2$")
	#plt.text(0.6,0.9,r"$m_\chi = 20$ GeV",transform=ax.transAxes)
plt.axhline(y=0.12,ls='dashed',c='k',label=r"$\Omega h^2 = 0.12$")
plt.legend()
plt.show()



### plot contour plot of relic abundance

###make plots pretty
plt.rc("text",usetex=True)
plt.rcParams.update({'mathtext.default':  'regular' })

plt.figure()
ax = plt.gca()
leveli = np.floor(np.log10(np.min(np.where(Omega > 0, Omega, np.inf))))
levelf = np.ceil(np.log10(np.max(Omega)))	
levels=np.power(10,np.arange(leveli,levelf+1,0.1))
contour = plt.contourf(np.log10(mdml),np.log10(yl),np.log10(Omega),levels=np.log10(levels))
contour_line = plt.contour(np.log10(mdml),np.log10(yl),np.log10(Omega),[np.log10(0.12)],colors='black',linestyles="solid") #contour at the correct abundance

#colour bar
cbar=plt.colorbar(contour)
cbarticks = np.arange(leveli,levelf+2,1)
cbarlabels = ['$10^{%d}$'%(v) for v in np.arange(leveli,levelf+2,1)]
cbar.set_ticks(cbarticks)
cbar.set_ticklabels(cbarlabels)
cbar.ax.set_ylabel(r'$\Omega h^2$',fontsize='12')
cbar.ax.tick_params(labelsize=12)

#axes
xticks = np.arange(np.ceil(np.log10(mdml[0])),np.log10(mdml[-1])+1,1)
yticks = np.arange(np.ceil(np.log10(yl[0])),np.log10(yl[-1]),2)
xlabels = ['$10^{%d}$'%(v) for v in xticks]
ylabels = ['$10^{%d}$'%(v) for v in yticks]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xlabels)
ax.set_yticklabels(ylabels)
plt.xlabel(r"$m_\chi\:(GeV)$")
plt.ylabel(r"$y = \epsilon^2 \alpha_D \left( \frac{m_\chi}{m_{A^\prime}}\right)^4$")
#plt.ylabel(r"$g_{SM}$")
plt.title(r"Vector DM with $A^\prime$ Mediator, $\alpha_D=%.1f$"%alphaD)

### add elastic scalar, pseudo-dirac
mdm_scalar,y_scalar = np.loadtxt("RelicAbundance/digitized_datasets/scalarelastic.txt",unpack=True,delimiter=",")
mdm_dirac,y_dirac = np.loadtxt("RelicAbundance/digitized_datasets/pseudodirac.txt",unpack=True,delimiter=",")
mdm_LDMX,y_LDMX = np.loadtxt("RelicAbundance/digitized_datasets/LDMX_projection.txt",unpack=True,delimiter=",")
plt.plot(np.log10(mdm_scalar/1000),np.log10(y_scalar),label=r"Scalar Elastic DM",color='tab:purple',ls="dashed")
plt.plot(np.log10(mdm_dirac/1000),np.log10(y_dirac),label=r"Pseudo-Dirac DM",color='black',ls="dotted")
plt.plot(np.log10(mdm_LDMX/1000),np.log10(y_LDMX),label="LDMX Projection",c='r')


contour_legend_lines = [Line2D([0], [0], color='black'),Line2D([0], [0], color='black',ls='dashed'),Line2D([0], [0], color='black',ls='dotted'),Line2D([0], [0], color='red')]
contour_legend_labels = ["Vector DM","Scalar Elastic DM", "Pseudo-Dirac DM", "LDMX Projection"]

plt.legend(contour_legend_lines,contour_legend_labels,loc="lower right")
# plt.ylim(-16,-6.8)
# plt.xlim(-3,np.log10(max(mdm_dirac/1000)))
#plt.savefig("RelicAbundance/alphaD0.5_g5_contour.pdf")
plt.show()






# ### plot <sigma v> (xf)
# plt.figure()
# for k in range(len(yl)):
# 	if k % 5 == 0.0:
# 		plt.loglog(mdml,sigmav_xf[k],label=r"$y = %.1e$"%yl[k])
# plt.legend()
# plt.xlabel(r"$m_\chi$ (GeV)")
# plt.ylabel(r"$<\sigma v>(x_f)$")

# # plt.figure()
# # for k in range(len(yl)):
# # 	if k % 5 == 0.0:
# # 		plt.semilogx(mdml,xf[k],label=r"$y = %.1e$"%yl[k])
# # plt.legend()
# # plt.xlabel(r"$m_\chi$ (GeV)")
# # plt.ylabel(r"$x_f = m_\chi / T_{FO}$")
# plt.show()


