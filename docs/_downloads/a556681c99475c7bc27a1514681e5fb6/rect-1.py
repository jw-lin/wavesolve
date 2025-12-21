from wavesolve.waveguide import Rectangle,Waveguide

airclad = Rectangle(1.,"air")
airclad.make_points(-12,12,-16,8)

strip1 = Rectangle(1.45,"strip1")
strip1.make_points(-5.5,-0.5,0,2)

strip2 = Rectangle(1.45,"strip2")
strip2.make_points(0.5,5.5,0,2)

corelayer = Rectangle(1.77,"corelayer")
corelayer.make_points(-10,10,-1,0)

clad = Rectangle(1.45,"clad")
clad.make_points(-10,10,-10,-1)

strip1.mesh_size = 1.0
strip2.mesh_size = 1.0
clad.mesh_size = 1.0
corelayer.mesh_size = 0.25

wvg = Waveguide([airclad,[strip1,strip2,corelayer,clad]] )
m = wvg.make_mesh(order=1)

wvg.plot_mesh(m,plot_points=False)