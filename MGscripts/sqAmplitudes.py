import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spyint

mf = 0.00511
h3 = b5 = g4=1.0

def gamma(mf,mAp):
	return (h3**2)/3 * (4*mAp**2+8*mf**2)

def scalar_amp_sq(th,s,mdm,mAp):
	return 2*h3**2*g4**2*(s-mdm**2)*(s-(s-4*mf**2)*np.cos(th)**2)/(s**2-mAp**4-2*mAp**2*s-mAp**2*gamma(mf,mAp)**2) * np.sin(th)

def scalar_dcs(s,mdm,mAp): # cross section (integrated squared amplitude * sin(th) over theta)
 	return spyint.quad(scalar_amp_sq,0,np.pi,args=(s,mdm,mAp))[0]

def vector_amp_sq(th,s,mdm,mAp):
	return b5**2*h3**2/(18*mdm**4)*(s-4*mdm**2)*(s**2-4*s*mdm**2+12*mdm**4)*(s-(s-4*mf**2)*np.cos(th)**2)/(s**2-mAp**4-2*mAp**2*s-mAp**2*gamma(mf,mAp)**2) * np.sin(th)

def vector_dcs(s,mdm,mAp):
 	return spyint.quad(vector_amp_sq,0,np.pi,args=(s,mdm,mAp))[0]

mdml = [1e-3,1e-2,1e-1]
colors = ['b','g','m']
j=0
for mdm in mdml:
	mAp = mdm*3
	sl = np.logspace(np.log10(4.1*mAp),np.log10(1e8*mAp),100)
	scalar = []
	vector = []
	for s in sl:
		scalar.append(scalar_dcs(s,mdm,mAp))
		vector.append(vector_dcs(s,mdm,mAp))
	plt.loglog(sl,scalar,label="scalar:mdm=%.2e"%mdm,c=colors[j],ls='dashed')
	plt.loglog(sl,vector,label="vector:mdm=%.2e"%mdm,c=colors[j])
	j+=1
plt.legend()
plt.show()