offsetsx = [1,np.cos(2*np.pi/3),np.cos(4*np.pi/3)]
offsetsy = [0,np.sin(2*np.pi/3),np.sin(4*np.pi/3)]
offsets = np.array([offsetsx,offsetsy]).T * core_offset

# make overlapping circular claddings
core_circs = []
for o in offsets:
    circ = waveguide.Circle(ncore,"") # label doesnt matter since we will combine these
    circ.make_points(rcore,res,o) # make points corresponding to radius rcore, resolution res, centerpoint o
    core_circs.append(circ)

# combine the Circles using Prim2DUnion, a type of Prim2D
core = waveguide.Prim2DUnion(core_circs,"core")

clad = waveguide.Circle(nclad,"clad")
clad.make_points(rclad,int(res/2)) # lower resolution for the outer boundary

PL3 = waveguide.Waveguide([clad,core])
PL3.plot_mesh()