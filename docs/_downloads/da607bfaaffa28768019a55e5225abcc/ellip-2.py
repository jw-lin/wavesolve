from wavesolve.FEsolver import count_modes,plot_vector_field,solve_waveguide_vec,get_eff_index
import matplotlib.pyplot as plt, numpy as np

wl = 1.55 # wavelength
IOR_dict = ellipfiber.assign_IOR()
ws,vs = solve_waveguide_vec(mesh,wl,IOR_dict,Nmax=4)
ne = get_eff_index(wl,ws) # convert all eigenvals into effective indices

# plotting
fig,axs = plt.subplots(2,2,sharex=True,sharey=True,figsize=(8,8))
for i,ax in enumerate(axs.flatten()):
    ax.set_aspect('equal')
    ax.set_title(r"$n_{\rm eff}=$"+str(round(np.real(ne[i]),5)))
    plot_vector_field(mesh,vs[i],show_mesh=False,ax=ax,bounds=(-5,5,-5,5))
    ellipfiber.plot_boundaries(ax) # show the core-cladding boundary
plt.subplots_adjust(hspace=0)
plt.show()