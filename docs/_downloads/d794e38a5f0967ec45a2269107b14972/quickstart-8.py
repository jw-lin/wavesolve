from wavesolve.FEsolver import plot_vector_field

fig,axs = plt.subplots(2,3,sharex=True,sharey=True,figsize=(12,8))
for i,ax in enumerate(axs.flatten()):
    plot_vector_field(mesh_ord1,eigvecs_ord1[i],ax=ax,bounds=(-8,8,-8,8))
plt.show()