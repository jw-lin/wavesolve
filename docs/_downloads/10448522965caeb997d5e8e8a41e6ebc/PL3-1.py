from wavesolve import waveguide
import numpy as np
import matplotlib.pyplot as plt

## params
ncore = 1.444+5e-3 # core index
nclad = 1.444 # cladding index

rcore = 4. # radius of circles composing flower-shaped core
rclad = 18. # cladding radius (outer simulation boundary)
res = 64 # resolution for boundaries

core_offset = 2. # radial offset of circles forming the core