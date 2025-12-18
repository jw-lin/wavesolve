from wavesolve.FEsolver import solve_waveguide

wl = 1.55 # um
eigvals, eigvecs = solve_waveguide(mesh,wl,IOR_dict)