from wavesolve.FEsolver import plot_scalar_field
import matplotlib.pyplot as plt

fig,axs = plt.subplots(1,3,sharey=True,figsize=(12,4))
for i,ax in enumerate(axs):
    ax.set_aspect('equal')
    plot_scalar_field(mesh,eigvecs[i],ax=ax,bounds=(-8,8,-8,8))
plt.show()