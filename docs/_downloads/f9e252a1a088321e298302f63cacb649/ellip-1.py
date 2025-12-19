from wavesolve import waveguide

## params
ncore = 1.512
nclad = 1.474
acore = 5.5/2 # core semimajor axis
bcore = 2.2/2 # core semiminor axis
rclad = 12. # outer boundary radius
core_res = 64 # resolution for core boundary
clad_res = 32 # resolution for outer boundary

## make the fiber
core = waveguide.Ellipse(ncore,"core")
core.make_points(acore,bcore,core_res)
clad = waveguide.Circle(nclad,"cladding")
clad.make_points(rclad,clad_res)
ellipfiber = waveguide.Waveguide([clad,core])

# make order 1 mesh
mesh = ellipfiber.make_mesh(order=1)

# only plot mesh lines
ellipfiber.plot_mesh(mesh,plot_points=False)