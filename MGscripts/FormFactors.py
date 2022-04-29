import matplotlib.pyplot as plt
import numpy as np

Me = 0.000511
mp = 0.938
mup = 2.79
Z = 74
Znuc=Z
A = 184
Anuc = A
scale_factor = 1.
inelastic1 = 1.9276/(scale_factor**2)
inelastic2 = 1.40845/(scale_factor**2)
aval = 111.0/(0.0005111*Znuc**(1./3.)) /scale_factor
apval = 773.0/(0.0005111*Znuc**(2./3.))/scale_factor
dval = 0.164/Anuc**(2./3.)*scale_factor**2


def G2_in(t):
	a = apval
	return (a**2*t/(1+a**2*t))**2 * ((1+t/(4*mp**2)*(mup**2-1))/(1+t/0.71)**4)**2 * Z

def G2_in_MG(t):
	return (apval**2*t/ (1+apval**2*t) )**2 * ((1+inelastic1*t)/ (1+inelastic2*t)**4 )**2 * Z

def G2_el(t):
	a = aval
	d = dval
	return (a**2 * t / (1+a**2*t))**2 * (1/(1+t/d))**2*Z**2

def G2_el_MG(t):
	return (aval **2 * t / (1+aval**2*t))**2 * (1 / (1+t/dval))**2 * Z**2



t = np.logspace(-5,5,100)


plt.loglog(t,G2_in(t),label="G2-inelastic",c='m')
plt.loglog(t,G2_in_MG(t),c='m',ls='--')
plt.loglog(t,G2_el(t),label="G2-elastic",c='g')
plt.loglog(t,G2_el_MG(t),c='g',ls='--')
plt.legend()
plt.show()
