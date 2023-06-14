""" use finite element method with quadratic triangular elements to solve for modes in the SCALAR approximation. 
    if this works, perhaps vector forthcoming. also, maybe fe-bpm?
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
#from numpy.linalg import eigh
from mesher import plot_mesh,lantern_mesh,circ_points,fiber_mesh,lantern_mesh_displaced_circles,lantern_mesh_3PL
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix

def IOR_fiber(r,n,n0):
    def inner(x,y):
        if x*x+y*y <= r*r:
            return n
        return n0
    return inner

def construct_AB(mesh,IOR_dict,k,poke_index = None):
    from shape_funcs import compute_dNdN, compute_NN
    
    points = mesh.points
    materials = mesh.cell_sets.keys()

    N = len(points)

    A = np.zeros((N,N))
    B = np.zeros((N,N))

    for material in materials:
        tris = mesh.cells[1].data[tuple(mesh.cell_sets[material])][0,:,0,:]
        for tri in tris:
            tri_points = points[tri]
            NN = compute_NN(tri_points)
            dNdN = compute_dNdN(tri_points)
            ix = np.ix_(tri,tri)
            A[ix] += (k**2*IOR_dict[material]**2) * NN - dNdN
            B[ix] += NN

    # now poke if necessary
    if poke_index is not None:
        for j,tri in enumerate(mesh.cells[1].data):
            if j == poke_index:
                tri_points = mesh.points[tri]
                NN = compute_NN(tri_points)
                ix = np.ix_(tri,tri)
                A[ix] += 0.1 * NN

    return A,B

def construct_AB_expl_IOR(mesh,IOR_arr,k):
    from shape_funcs import compute_dNdN, compute_NN
    
    points = mesh.points
    materials = mesh.cell_sets.keys()

    N = len(points)

    A = np.zeros((N,N))
    B = np.zeros((N,N))

    for tri,IOR in zip(mesh.cells[1].data,IOR_arr):
        tri_points = points[tri]
        NN = compute_NN(tri_points)
        dNdN = compute_dNdN(tri_points)
        ix = np.ix_(tri,tri)
        A[ix] += (k**2*IOR**2) * NN - dNdN
        B[ix] += NN

    return A,B

# turns out... this isn't actually faster. lol. completely dominated by solving system, not generating matrix
def construct_AB_fast(mesh,IOR_dict,k):
    from shape_funcs import compute_dNdN_precomp, compute_J_and_mat, compute_NN_precomp
    
    points = mesh.points
    verts = mesh.cells[1].data
    materials = mesh.cell_sets.keys()

    N = len(points)

    A = np.zeros((N,N))
    B = np.zeros((N,N))

    IORs = np.zeros(N)
    
    for material in materials:
        IORs[tuple(mesh.cell_sets[material])] = IOR_dict[material]

    _Js,mat = compute_J_and_mat(points,verts)

    for _J,row,IOR,idxs in zip(_Js,mat,IORs,verts):
        NN = compute_NN_precomp(_J)
        dNdN = compute_dNdN_precomp(row)
        ix = np.ix_(idxs,idxs)
        A[ix] += k**2*IOR**2 * NN - dNdN
        B[ix] += NN

    return A,B

def solve(A,B,mesh,k,IOR_dict,plot=False):
    w,v = eigh(A,B)

    IORs = [ior[1] for ior in IOR_dict.items()]
    nmin,nmax = min(IORs) , max(IORs)
    mode_count = 0
    
    for _w,_v in zip(w[::-1],v.T[::-1]):
        if _w<0:
            continue
        ne = np.sqrt(_w/k**2)
        if plot:
            if not (nmin <= ne <= nmax):
                print("warning: spurious mode! ")
            
            print("effective index: ",np.sqrt(_w/k**2))
            plot_eigenvector(mesh,_v)
        if (nmin <= ne <= nmax):
            mode_count+=1
        else:
            break
    
    # normalization
    # v /= np.sqrt(np.sum(np.power(v,2),axis=0))[None,:]

    return w[::-1],v.T[::-1],mode_count

def solve_sparse(A,B,mesh,k,IOR_dict,plot=False,num_modes=6):
    w,v = eigsh(A,M=B,k=num_modes,which="LA")
    IORs = [ior[1] for ior in IOR_dict.items()]
    nmin,nmax = min(IORs) , max(IORs)
    mode_count = 0
    
    for _w,_v in zip(w[::-1],v.T[::-1]):
        if _w<0:
            continue
        ne = np.sqrt(_w/k**2)
        if plot:
            if not (nmin <= ne <= nmax):
                print("warning: spurious mode! ")
            
            print("effective index: ",np.sqrt(_w/k**2))
            plot_eigenvector(mesh,_v)
        if (nmin <= ne <= nmax):
            mode_count+=1
        else:
            break
    
    # normalization
    # v /= np.sqrt(np.sum(np.power(v,2),axis=0))[None,:]

    return w[::-1],v.T[::-1],mode_count

def plot_eigenvector(mesh,v,plot_mesh = False,plot_circle=False):
    points = mesh.points
    fig,ax = plt.subplots(figsize=(5,5))
    plt.axis('equal')
    plt.tricontourf(points[:,0],points[:,1],v,levels=60)
    plt.colorbar()
    if plot_mesh:
        plot_mesh(mesh,show=False,ax=ax)
    if plot_circle:
        circle = plt.Circle((0,0),mesh.cell_data['radius'],ec='white',fc='None',lw=2)
        ax.add_patch(circle)
    plt.show()

def compute_diff(tri_idx,mesh,_pinv):
    from shape_funcs import compute_NN
    point_idxs = mesh.cells[1].data[tri_idx]
    points = mesh.points[point_idxs]
    N = len(mesh.points)
    ix = np.ix_(point_idxs,point_idxs)
    B_tri = compute_NN(points)

    return np.dot(_pinv[ix],-B_tri),ix,point_idxs

def compute_IOR_arr(mesh,IOR_dict):
    IORs = np.zeros(len(mesh.cells[1].data))
    materials = mesh.cell_sets.keys()
    for material in materials:
        IORs[tuple(mesh.cell_sets[material])] = IOR_dict[material]
    
    return IORs

def optimize_for_mode_structure(mesh,IOR_dict,k,target_field,iterations = 1):
    '''find the refractive index profile such that the fundamental mode most closely matches the target field'''

    # first, solve the base mesh
    A,B = construct_AB(mesh,IOR_dict,k)
    w,v = solve(A,B,mesh,k,IOR_dict,plot=False)

    for _it in range(iterations):

        # next, compute the matrix of derivatives
        _pinv = np.linalg.pinv(A - w[0]*B)

        N_tris = len(mesh.cells[1].data)
        N_points = A.shape[0]

        diff_mat = np.zeros((N_tris,N_points))

        for i in range(N_tris):
            _diff,ix,point_idxs = compute_diff(i,mesh,_pinv)
            diff_field = np.dot(_diff,v[0][point_idxs])
            diff_field /= np.sqrt(np.sum(np.power(diff_field,2)))
            diff_mat[i,point_idxs] = diff_field

        # next, compute the difference between v[0] and the target
        diff = v[0] - target_field

        # compute the overlaps between the derivatives and the difference field
        coeffs = np.dot(diff_mat,diff)

        # subtract off the coeffs from refractive index profile
        IOR0 = compute_IOR_arr(mesh,IOR_dict)
        IOR = np.sqrt(np.power(IOR0,2)*k**2-coeffs)/k

        # check the new results
        A,B = construct_AB_expl_IOR(mesh,IOR,k)
        w,v = solve(A,B,mesh,k,IOR_dict)

    # plot the result
    #plot_eigenvector(IOR)
    from mesher import plot_mesh_expl_IOR
    plot_mesh_expl_IOR(mesh,IOR)

    plot_eigenvector(mesh,v[0])

    xcs = []
    ycs = []
    for tri in mesh.cells[1].data:
        tri_points = mesh.points[tri]
        xm = np.mean(tri_points[:,0])
        ym = np.mean(tri_points[:,1])
        xcs.append(xm)
        ycs.append(ym)

    plt.tricontourf(xcs,ycs,IOR,levels=40)
    plt.colorbar()
    plt.show()

"""
ncore = 1.444 + 0.01036
nclad = 1.444
njack = 1.444 - 5.5e-3

k = 2*np.pi/1.55
res = 40

r_clad = 10
w = 4*r_clad
r_core = 2.2/4

### core shift analysis

cores = [(0,0)] + circ_points(2*r_clad/3,5)
IOR_dict = {"core0":ncore,"core1":ncore,"core2":ncore,"core3":ncore,"core4":nclad,"core5":nclad,"core6":nclad,"cladding":nclad,"jacket":njack}

IOR_dict2 = IOR_dict = {"core0":ncore,"core1":ncore,"core2":nclad,"core3":nclad,"core4":ncore,"core5":ncore,"core6":nclad,"cladding":nclad,"jacket":njack}

mesh = lantern_mesh_displaced_circles(w/2,r_clad,cores,r_core,30,ds=0.2)
#mesh = lantern_mesh(w/2,r_clad,cores,r_core,res,petals=5,petal_amp=0.2)

print("mesh and refractive index distribution")
plot_mesh(mesh,IOR_dict=IOR_dict)

A,B = construct_AB(mesh,IOR_dict,k)
w,v = solve(A,B,mesh,k,IOR_dict,plot=False)


_A,_B = construct_AB(mesh,IOR_dict2,k)
_w,_v = solve(_A,_B,mesh,k,IOR_dict2,plot=False)

for vvec, _vvec in zip(v,_v):
    plot_eigenvector(mesh,_vvec-vvec)
"""

# target mode, v[5]

#optimize_for_mode_structure(mesh,IOR_dict,k,v[1],iterations=1)


def get_num_modes(wl,r):
    k = 2*np.pi/wl
    ncore = 1.444
    nclad = 1.444-5.5e-3
    m = lantern_mesh_3PL(r,16)
    IOR_dict = {"jacket":nclad,"cladding":ncore}
    A,B = construct_AB(m,IOR_dict,k)
    _A = csr_matrix(A)
    _B = csr_matrix(B)

    w,v,n = solve_sparse(_A,_B,m,k,IOR_dict,plot=False)
    return n


wl_arr = np.linspace(1,1.8,20)
r_arr = np.linspace(3.3-0.5,5.3-0.5,20)

out = np.zeros((20,20))
ncore = 1.444
nclad = 1.444-5.5e-3
IOR_dict = {"jacket":nclad,"cladding":ncore}
itot = 0
for j in range(20):
    m = lantern_mesh_3PL(r_arr[j],16)
    for i in range(20):
        print("iteration "+str(itot))
        wl = wl_arr[i]
        k = 2*np.pi/wl
        A,B = construct_AB(m,IOR_dict,k)
        #_A = csr_matrix(A)
        #_B = csr_matrix(B)
        w,v,n = solve(A,B,m,k,IOR_dict,plot=False)
        out[i,j] = n
        itot+=1

np.save("countmap",out)
plt.imshow(out.T,origin='lower',extent=(1,1.8,2.8*2/np.sqrt(3),4.8*2/np.sqrt(3)),aspect=4/np.sqrt(3)/0.8)
plt.xlabel("wavelength")
plt.ylabel("core radius")
plt.show()
