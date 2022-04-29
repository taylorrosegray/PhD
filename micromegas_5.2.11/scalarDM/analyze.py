import matplotlib.pyplot as plt
import numpy as np

mAp = 3.0
T,sigmav = np.loadtxt("sigmav_mAp%.1e.dat"%mAp,unpack=True)

plt.loglog(T,sigmav)
plt.xlabel("T (GeV)")
plt.ylabel(r"$<\sigma v>$")
plt.show()