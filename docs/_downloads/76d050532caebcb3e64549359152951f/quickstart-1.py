from wavesolve.waveguide import CircularFiber

#params
rcore = 5               # core radius, assumed um
rclad = 15              # outer radius of simulation boundary, um
ncore = 1.444+8.8e-3    # core refractive index
nclad = 1.444           # cladding refractive index
core_res = 64           # number of line segments used to refine the core-cladding boundary
clad_res = 64           # number of line segments for outer clad boundary

# make waveguide
fiber = CircularFiber(rcore,rclad,ncore,nclad,core_res,clad_res,core_mesh_size=0.5)

# make a plot mesh
mesh = fiber.make_mesh()
fiber.plot_mesh(mesh=mesh)