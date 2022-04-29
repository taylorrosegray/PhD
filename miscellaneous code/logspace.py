import numpy as np
import matplotlib.pyplot as plt

def f(i,a,b,k):
	return 10**(a+i*(b-a)/(k))

a,b = -3.0,1.0
k = 10
i = np.linspace(0,10,k)

plt.scatter(i,f(i,a,b,k),label="my function")
plt.plot(i,np.logspace(a,b,k),label="numpy")
plt.legend()
plt.show()

