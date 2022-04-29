import numpy as np
import matplotlib.pyplot as plt
import os

def eps_f(y,alphad,mdm,mAp):
	return np.sqrt(y/alphad * (mAp/mdm)**4)

### input parameters
alphad_value = 0.5
gD_value = np.sqrt(4*np.pi*alphad_value)

epsilon = 0.01
epsilon_nu = 0.0
epsilon_e = 0.0
epsilon_mu = 0.0
epsilon_tau = 0.0
epsilon_d = 0.0
epsilon_u = 0.01

epsilon_ul = [0.0,0.01]
epsilon_el = [0.01,0.0]
lsl = ['-','--']

y = 1e-7
input_file = "data_python.par"

### loop over
#mdml = [1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,1e-2]
#mdml=np.logspace(-3,0,10)
mdml= [100.0]

cmap = plt.get_cmap('gnuplot')
colors = [cmap(i) for i in np.linspace(0, 1, len(mdml))]

Xfl = []

#for i in range(len(epsilon_el)):
mdmi = 0
sigmav = []
for mdm in mdml:
	
	mAp = mdm*3
	#epsilon = eps_f(y,alphad_value,mdm,mAp)

	print("mdm = %.2e" %mdm)

	### write input file
	Mdm    = "mdm    %.1E" % mdm
	MAp    = "mAp    %.1E" % mAp
	eps_nu = "eps_nu %.3E" % epsilon
	eps_e  = "eps_e %.3E" % epsilon
	eps_mu = "eps_mu %.3E" % epsilon
	eps_tau= "eps_tau %.3E" % epsilon
	eps_d  = "eps_d %.3E" % epsilon
	eps_u  = "eps_u %.3E" % epsilon
	gD     = "gD     %.2E"  %gD_value

	input_parameters = [Mdm,MAp,eps_nu,eps_e,eps_mu,eps_tau,eps_d,eps_u,gD]
	np.savetxt(input_file,input_parameters, delimiter=" ", fmt="%s")

	### run main.c
	os.system("./mymain " + input_file)

	### read main output file
	T_file,sigmav_file = np.loadtxt("sigmav.dat",unpack=True)
	#sigmav.append(sigmav_file)
	#Xf,omega_file = np.loadtxt("Omega.dat",unpack=True)
	#Xfl.append(Xf)

	plt.loglog(T_file,sigmav_file,label=r"$m_{\chi} = %.1e GeV$"%mdm,c=colors[mdmi])
	plt.axvline(x=mdm,c=colors[mdmi],ls=':')
	#plt.axvline(x=mdm/Xf, ls='--',c=colors[mdmi])

	mdmi +=	1

#plt.text(1e4,0.01,r"$\epsilon = 10^{%1d}$"%np.log10(epsilon))
#plt.loglog(mdml,sigmav)
#plt.ylabel(r"$\sigma v(m_\chi)$")
#plt.legend()
#plt.title(r"$\chi \bar\chi \to u u$")
plt.title("sigmav")
plt.xlabel("T")
plt.show()


