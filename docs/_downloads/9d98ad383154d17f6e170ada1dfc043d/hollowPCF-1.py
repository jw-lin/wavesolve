from wavesolve.waveguide import PhotonicBandgapFiber

void_radius = 11.5/2    # hollow core radius
hole_radius = 2.        # radius for air holes in cladding
clad_radius = 15.       # outer boundary radius
nclad = 1.444           # cladding index
hole_separation = 4.05  # separation between air holes

hollow_PCF = PhotonicBandgapFiber(void_radius,hole_radius,clad_radius,nclad,hole_separation,20,32,hole_mesh_size=1.2,clad_mesh_size=2.)
m = hollow_PCF.make_mesh(order=1)
hollow_PCF.plot_mesh(m)