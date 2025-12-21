from wavesolve.waveguide import PhotonicBandgapFiber

void_radius = 11.5/2    # hollow core radius
hole_radius = 2.        # radius for air holes in cladding
clad_radius = 15.       # outer boundary radius
nclad = 1.444           # cladding index
hole_separation = 4.05  # separation between air holes
hole_res = 24           # air hole boundary resolution
clad_res = 32           # outer boundary (cladding) resolution

hollow_PCF = PhotonicBandgapFiber(void_radius,hole_radius,clad_radius,nclad,hole_separation,hole_res,clad_res,hole_mesh_size=1.0,clad_mesh_size=2.0)
m = hollow_PCF.make_mesh(order=1)
hollow_PCF.plot_mesh(m)