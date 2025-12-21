from wavesolve.FEsolver import get_eff_index,plot_vector_field,solve_waveguide_vec
import matplotlib.pyplot as plt, numpy as np

IOR_dict = hollow_PCF.assign_IOR()
wl = 1.65

ws,vs = solve_waveguide_vec(m,wl,IOR_dict,target_neff=1.0,Nmax=10)
ne = get_eff_index(wl,ws)

fig,axs = plt.subplots(2,3,sharey=True,figsize=(10,6))
for j,ax in enumerate(axs.T.flatten()):
    i = j+3
    plot_vector_field(m,vs[i],ax=ax,bounds=(-8,8,-8,8))
    hollow_PCF.plot_boundaries(ax)
    ax.set_title("mode "+str(i+1)+r": $n_{\rm eff}=$"+str(round(np.real(ne[i]),4)))

plt.show()