import matplotlib.pyplot as plt
from wavesolve.FEsolver import solve_waveguide, plot_scalar_field, count_modes

IOR_dict = PL3.assign_IOR()
mesh = PL3.make_mesh()
wl = 1.55

# not sure how many modes there will be so we'll use the
# optional Nmax arg to get the first 10 modes
_w,_v = solve_waveguide(mesh,wl,IOR_dict,Nmax=10)

# get number of modes
Nmodes = count_modes(_w,wl,IOR_dict)

fig,axs = plt.subplots(1,3,sharex=True,sharey=True,figsize=(12,4))
plot_scalar_field(mesh,_v[0],show_mesh=False,ax=axs[0])
plot_scalar_field(mesh,_v[Nmodes-1],show_mesh=False,ax=axs[1])
plot_scalar_field(mesh,_v[Nmodes],show_mesh=False,ax=axs[2])

for ax in axs:
    PL3.plot_boundaries(ax) # show the waveguide boundary

axs[0].set_title("mode 1")
axs[1].set_title("mode "+str(Nmodes))
axs[2].set_title("mode "+str(Nmodes+1))

plt.show()