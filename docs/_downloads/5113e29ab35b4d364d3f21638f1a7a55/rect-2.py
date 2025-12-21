from wavesolve.FEsolver import solve_waveguide_vec,plot_vector_field,get_eff_index
import matplotlib.pyplot as plt

wl = 1.
IOR_dict = wvg.assign_IOR()

w,v = solve_waveguide_vec(m,wl,IOR_dict,Nmax = 8)

ne = get_eff_index(wl,w)

fig,axs = plt.subplots(8,1,sharex=True,sharey=True,figsize=(8,12))
for i,ax in enumerate(axs.flatten()):
    ax.set_title("mode "+str(i+1)+r": $n_{\rm eff}=$"+str(round(np.real(ne[i]),5)))
    plot_vector_field(m,v[i],ax=axs.flatten()[i],bounds=(-10,10,-2,1))
    wvg.plot_boundaries(ax)
plt.show()