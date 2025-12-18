from wavesolve.FEsolver import count_modes
num_modes = count_modes(eigvals,wl,IOR_dict)
print("number of guided modes: ",num_modes)