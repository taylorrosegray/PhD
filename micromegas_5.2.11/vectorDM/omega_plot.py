import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

alphad = 0.5
#m_ratiol = [3.0,10.0] # mAp/mdm
m_ratio = 3.0
cases = ["b5","b6R","b6I","b7R","b7I"]

rd = 0.12
colors = ["tab:blue","tab:orange","tab:green","tab:cyan","tab:purple"]
contour_legend_lines = [Line2D([0], [0], color='black',ls='dashed'),Line2D([0], [0], color='black',ls='dotted'),Line2D([0], [0], color='red')]
contour_legend_labels = ["Scalar Elastic DM", "Pseudo-Dirac DM", "LDMX Projection"]

plt.figure()
ax = plt.gca()
r=0
#for m_ratio in m_ratiol:
for case in cases:
	run_name = "vectorDM_"+case+"_ratio"+str(m_ratio)+"_alphaD"+str(alphad)
	path = "runs/"+run_name+"/"

	#import data
	omegal = np.load(path+"omega.npy")
	yl = np.load(path+"y.npy")
	mdml = np.load(path+"mdm.npy")

	### plot
	CL  = plt.contour(np.log10(mdml),np.log10(yl),np.log10(omegal),[np.log10(rd)],colors=colors[r],linestyles='solid')
	contour_legend_lines.append(Line2D([0], [0], color=colors[r]))
	#contour_legend_labels.append(r"$m_{A^\prime}/m_\chi=%d$"%m_ratio)
	contour_legend_labels.append("Vector DM:"+case)
	r+=1

plt.rc("font",family="serif",serif="Computer Modern")
plt.rc("text",usetex=True)
plt.rcParams.update({'mathtext.default': 'regular'})

#add elastic scalar, pseudo-dirac
mdm_scalar,y_scalar = np.loadtxt("digitized_datasets/scalarelastic.txt",unpack=True,delimiter=",")
mdm_dirac,y_dirac = np.loadtxt("digitized_datasets/pseudodirac.txt",unpack=True,delimiter=",")
mdm_LDMX,y_LDMX = np.loadtxt("digitized_datasets/LDMX_projection.txt",unpack=True,delimiter=",")
plt.plot(np.log10(mdm_scalar/1000),np.log10(y_scalar),label=r"Scalar Elastic DM",color='black',ls="dashed")
plt.plot(np.log10(mdm_dirac/1000),np.log10(y_dirac),label=r"Pseudo-Dirac DM",color='black',ls="dotted")
plt.plot(np.log10(mdm_LDMX/1000),np.log10(y_LDMX),label="LDMX Projection",c='red')

xticks = np.arange(np.ceil(np.log10(mdml[0])),np.log10(mdml[-1]),1)
yticks = np.arange(np.ceil(np.log10(yl[0])),np.log10(yl[-1]),1)
xlabels = ['$10^{%d}$'%(v) for v in xticks]
ylabels = ['$10^{%d}$'%(v) for v in yticks]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xlabels)
ax.set_yticklabels(ylabels)
plt.xlabel(r"$m_\chi$ (GeV)")
plt.ylabel(r"$y = \epsilon^2 \alpha_D \left( \frac{m_\chi}{m_{A^\prime}}\right)^4$")
plt.ylim(-16,-7)
plt.xlim(-3,np.log10(max(mdml)))
plt.text(0.0,-11.5,r"$\alpha_D = %.1f$"%alphad)
plt.legend(contour_legend_lines,contour_legend_labels,loc="lower right")

#lt.axvline(x=np.log10(0.00255),c='k')
plt.axvline(x=np.log10(0.1),c='k',ls='dashed')

#plt.savefig("VectorDM_varycouplingtypes.pdf")	
plt.show()





