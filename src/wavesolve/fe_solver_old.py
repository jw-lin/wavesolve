""" use finite element method with quadratic triangular elements to solve for modes in the SCALAR approximation. 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh,eig
from scipy.sparse.linalg import eigsh,eigs,spsolve
from scipy.sparse import csr_matrix, lil_matrix
from wavesolve.shape_funcs import affine_transform, get_basis_funcs_affine,apply_affine_transform,evaluate_basis_funcs,get_linear_basis_funcs_affine,get_edge_linear_basis_funcs_affine
from wavesolve.mesher import construct_meshtree
from wavesolve.shape_funcs import *
from wavesolve.shape_funcs import compute_NN_dNdN_vec
from wavesolve.waveguide import plot_mesh
try:
    import pypardiso
except:
    pass

#region FEM matrices

def construct_AB_order2(mesh,IOR_dict,k,sparse=False,poke_index = None):
    """ construct the A and B matrices corresponding to the given waveguide geometry.
    Args:
    mesh: the waveguide mesh, produced by wavesolve.mesher or pygmsh
    IOR_dict: dictionary that assigns each named region in the mesh to a refractive index value
    k: free-space wavenumber of propagation

    Returns:
    A,B: matrices
    """
    
    points = mesh.points
    materials = mesh.cell_sets.keys()

    N = len(points)
        
    all_tris = [mesh.cells[1].data[tuple(mesh.cell_sets[material])][0,:,0,:] for material in materials]
    if sparse:
        all_tris_v = np.vstack(all_tris)
        row_indices = np.repeat(all_tris_v, 6, 0).flatten()
        col_indices = np.repeat(all_tris_v, 6, 1).flatten()
        A_data = np.zeros(len(row_indices))
        B_data = np.zeros(len(row_indices))
    else:
        A = np.zeros((N,N))
        B = np.zeros((N,N))
        
    mat_idx = 0
    for (tris, material) in zip(all_tris, materials):
        tris_points = points[tris]
        NNs, dNdNs = compute_NN_dNdN_vec(tris_points)
        A_places = ((k**2*IOR_dict[material]**2) * NNs - dNdNs).flatten()
        mat_N = len(A_places)
        A_data[mat_idx:mat_idx+mat_N] += A_places
        B_data[mat_idx:mat_idx+mat_N] += NNs.flatten()
        mat_idx += mat_N

    # now poke if necessary
    if poke_index is not None:
        for j,tri in enumerate(mesh.cells[1].data):
            if j == poke_index:
                tri_points = mesh.points[tri]
                NN = compute_NN(tri_points)
                ix = np.ix_(tri,tri)
                A[ix] += 0.1 * NN

    if sparse:
        A = csr_matrix((A_data, (row_indices, col_indices)))
        B = csr_matrix((B_data, (row_indices, col_indices)))
    return A,B

def construct_B_order2(mesh,sparse=False):
    """ construct only the B matrix ("mass matrix") corresponding to the given waveguide geometry. this is used for inner products.
    Args:
    mesh: the waveguide mesh, produced by wavesolve.mesher or pygmsh
    k: free-space wavenumber of propagation

    Returns:
    B: matrix
    """
    
    points = mesh.points
    tris = mesh.cells[1].data
    tris_v = np.vstack(tris)
    tris_points = points[tris]
    NNs, dNdNs = compute_NN_dNdN_vec(tris_points)
    if sparse:
        row_indices = np.repeat(tris_v, 6, 0).flatten()
        col_indices = np.repeat(tris_v, 6, 1).flatten()
        B = csr_matrix((NNs.flatten(), (row_indices, col_indices)))
    else:
        N = len(points)
        B = np.zeros((N,N))
        for tri in mesh.cells[1].data:
            tri_points = points[tri]
            NN = compute_NN(tri_points)
            ix = np.ix_(tri,tri)
            B[ix] += NN

    return B

construct_B = construct_B_order2

def construct_AB_order1(mesh,IOR_dict,k,sparse=False):
    points = mesh.points
    tris = mesh.cells[1].data 
    materials = mesh.cell_sets.keys()

    N = len(points)
    if not sparse:
        A = np.zeros((N,N))
        B = np.zeros((N,N))
    else:
        A = lil_matrix((N,N))
        B = lil_matrix((N,N))

    for material in materials:
        tris = mesh.cells[1].data[tuple(mesh.cell_sets[material])][0,:,0,:]

        for tri in tris:
            tri_points = points[tri]
            pc = precompute(tri_points,tri)

            NN = computeL_NN(precomp=pc)
            dNdN = computeL_dNdN(precomp=pc)

            ix = np.ix_(tri,tri)
            A[ix] += (k**2*IOR_dict[material]**2) * NN - dNdN
            B[ix] += NN
    return A,B

def construct_AB(mesh,IOR_dict,k,sparse=False,order=2):
    """ construct the generalized eigenvalue problem matrices for the SCALAR formulation of FEM
    ARGS
        mesh: the finite element mesh object
        IOR_dict: dictionary of refractive index values corresponding to mesh
        k: free space wavenumber
        sparse: whether to return the matrices as sparse or dense; sparse is typically a better option for large meshes with ~>1000 nodes
    """
    if order == 2:
        assert mesh.cells[1].data.shape[1] == 6, "must use order 2 mesh for order 2 solver."
        return construct_AB_order2(mesh,IOR_dict,k,sparse)
    elif order == 1:
        assert mesh.cells[1].data.shape[1] == 3, "must use order 1 mesh for order 1 solver"
        return construct_AB_order1(mesh,IOR_dict,k,sparse)
    else:
        raise NotImplementedError

def construct_AB_vec(mesh,IOR_dict,k,sparse=False):
    """construct generalized eigenvalue problem matrices for VECTOR formulation of FEM
    ARGS
        mesh: the finite element mesh object
        IOR_dict: dictionary of refractive index values corresponding to mesh
        k: free space wavenumber
        sparse: whether to return the matrices as sparse or dense; sparse is typically a better option for large meshes with ~>1000 nodes
    """
    points = mesh.points
    tris = mesh.cells[1].data 
    materials = mesh.cell_sets.keys()
    edges = mesh.cells[0].data # for this to work, need to update mesh with get_unique_edges()
    Ntt = len(edges)
    Nzz = len(points)
    N = Ntt+Nzz
    if not sparse:
        Att = np.zeros((Ntt,Ntt))
        A = np.zeros((N,N))
        B = np.zeros((N,N))
        Bzz = np.zeros((Nzz,Nzz))
        Btz = np.zeros((Ntt,Nzz))
        Btt = np.zeros((Ntt,Ntt))
    else:
        Att = lil_matrix((Ntt,Ntt))
        A = lil_matrix((N,N))
        B = lil_matrix((N,N))
        Bzz = lil_matrix((Nzz,Nzz))
        Btz = lil_matrix((Ntt,Nzz))
        Btt = lil_matrix((Ntt,Ntt))

    for material in materials:
        tris = mesh.cells[1].data[tuple(mesh.cell_sets[material])][0,:,0,:]
        edge_indices = mesh.edge_indices[tuple(mesh.cell_sets[material])][0,:,0,:]
        _k2 = (k**2*IOR_dict[material]**2)
        for tri,idx in zip(tris,edge_indices):
            tri_points = points[tri]
            pc = precompute(tri_points,tri)

            NeNe = computeL_Ne_Ne(precomp=pc)
            NN = computeL_NN(precomp=pc)
            dNdN = computeL_dNdN(precomp=pc)
            NedN = computeL_Ne_dN(precomp=pc)
            cdNcdN = computeL_curlNe_curlNe(precomp=pc)

            ixtt = np.ix_(idx,idx)
            ixtz = np.ix_(idx,tri)
            ixzz = np.ix_(tri,tri)

            Att[ixtt] += _k2 * NeNe - cdNcdN
            Btt[ixtt] += NeNe
            Btz[ixtz] += NedN
            Bzz[ixzz] += dNdN - _k2*NN
    
    _ixtt = np.ix_(range(Ntt),range(Ntt))
    _ixtz = np.ix_(range(Ntt),range(Ntt,Ntt+Nzz))
    _ixzt = np.ix_(range(Ntt,Ntt+Nzz),range(Ntt))
    _ixzz = np.ix_(range(Ntt,Ntt+Nzz),range(Ntt,Ntt+Nzz))

    A[_ixtt] += Att
    B[_ixtt] += Btt
    B[_ixtz] += Btz
    B[_ixzt] += Btz.transpose()
    B[_ixzz] += Bzz

    if sparse:
        return A.tocsc(),B.tocsc()
    return A,B

#endregion

#region FEM solving

def solve(A,B,mesh,k,IOR_dict,plot=False):
    """ Given the A,B matrices, solve the general eigenvalue problem A v = w B v
        where v are the eigenvectors and w the eigenvalues.
        Args:
        A: A matrix of eigenvalue problem
        B: B matrix of eigenvalue
        mesh: the waveguide mesh, produced by wavesolve.mesher or pygmsh
        k: free-space wavenumber of propagation
        IOR_dict: dictionary that assigns each named region in the mesh to a refractive index value
        plot: set True to plot eigenvectors in descending order of eigenvalue

        Returns: 
        w: list of eigenvalues in descending order
        v: list of eigenvectors
        N: number of non-spurious eigenvectors (guided modes) found.
    """
    w,v = eigh(A,B,overwrite_a=True,overwrite_b=True)

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
    """An extension of solve() to A and B matrices in CSR format."""
    
    est_eigval = np.power(k*max(IOR_dict.values()),2)
    w,v = eigsh(A.tocsr(),M=B.tocsr(),k=num_modes,which="SA",sigma=est_eigval)

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

def solve_waveguide(mesh,wl,IOR_dict,plot=False,ignore_warning=False,sparse=True,Nmax=10,order=2,target_neff=None):
    """ given a mesh, propagation wavelength, and refractive index dictionary, solve for the SCALAR modes. 
        this has the same functionality as running construct_AB() and solve(). 
    
    ARGS: 
        mesh: mesh object corresponding to waveguide geometry
        wl: wavelength, defined in the same units as mesh point positions
        IOR_dict: a dictionary assigning different named regions of the mesh different refractive index values
        plot: set True to view eigenmodes
        ignore_warning: bypass the warning raised when the mesh becomes too large to solve safely with scipy.linalg.eigh()
        sparse: set True to use a sparse solver, which is can handle larger meshes but is slower
        Nmax: return only the <Nmax> largest eigenvalue/eigenvector pairs
        order: the order of the triangular finite elements. can be 1 (linear) or 2 (quadratic) ; default 2
        target_neff: search for modes with indices close to but below this value. if None, target_neff is set to the maximum index in the guide.
    RETURNS:
        w: array of eigenvalues, descending order
        v: array of corresponding eigenvectors (waveguide modes)
        N: number non-spurious (i.e. propagating) waveguide modes
    """
    
    k = 2*np.pi/wl
    if target_neff is None:
        est_eigval = np.power(k*max(IOR_dict.values()),2)
    else:
        est_eigval = np.power(k*target_neff,2)
    
    A,B = construct_AB(mesh,IOR_dict,k,sparse=sparse,order=order)
    N = A.shape[0]

    if A.shape[0]>2000 and not ignore_warning and not sparse:
        raise Exception("A and B matrices are larger than 2000 x 2000 - this may make your system unstable. consider setting sparse=True")
    if not sparse:
        w,v = eigh(A,B,subset_by_index=[N-Nmax,N-1],overwrite_a=True,overwrite_b=True)
    else:
        _A = A.tocsr()
        _B = B.tocsr()
        w,v = eigsh(_A,M=_B,k=Nmax,which="SA",sigma=est_eigval)

    IORs = [ior[1] for ior in IOR_dict.items()]
    nmin,nmax = min(IORs) , max(IORs)
    mode_count = 0
    
    for _w,_v in zip(w[::-1],v.T[::-1]):
        if _w<0:
            continue
        ne = np.sqrt(_w/k**2)
        if plot:
            if not (nmin <= ne <= nmax):
                print("warning: spurious mode! stopping plotting ... ")
            print("effective index: ",get_eff_index(wl,_w))
            plot_scalar_mode(mesh,_v)
        if (nmin <= ne <= nmax):
            mode_count+=1
        else:
            break

    return w[::-1],v.T[::-1],mode_count

def solve_waveguide_vec(mesh,wl,IOR_dict,plot=False,ignore_warning=False,sparse=True,Nmax=10,target_neff=None,sparse_solve_mode='transform',verbose=True):
    """ given a mesh, propagation wavelength, and refractive index dictionary, solve for VECTOR modes, using linear triangles (order 1).
    
    ARGS: 
        mesh: mesh object corresponding to waveguide geometry
        wl: wavelength, defined in the same units as mesh point positions
        IOR_dict: a dictionary assigning different named regions of the mesh different refractive index values
        plot: set True to view eigenmodes
        ignore_warning: bypass the warning raised when the mesh becomes too large to solve safely with scipy.linalg.eigh()
        sparse: set True to use a sparse solver, which is can handle larger meshes but is slower
        Nmax: return only the <Nmax> largest eigenvalue/eigenvector pairs
        target_neff: search for modes with indices close to but below this value. if None, target_neff is set to the maximum index in the guide.
        sparse_solve_mode: mode to solve the generalized eigenvalue problem. default is 'transform' -> use spsolve to convert to ordinary eigenvalue problem. 
        verbose: set False to block all printouts
    RETURNS:
        w: array of eigenvalues, descending order
        v: array of corresponding eigenvectors (waveguide modes)
        N: number non-spurious (i.e. propagating) waveguide modes
    """

    assert mesh.cells[1].data.shape[1] == 3, "must use order 1 mesh for vectorial solver"
    assert sparse_solve_mode in ["straight","transform","pardiso"], "sparse solve mode must be `straight` (plug into eigenproblem into eigsh) or `transform` (use spsolve to convert into an ordinary eigenproblen, then use eigs)"

    k = 2*np.pi/wl

    if target_neff is None:
        est_eigval = np.power(k*max(IOR_dict.values()),2)
    else:
        est_eigval = np.power(k*target_neff,2)

    if verbose:
        print("building matrices ... ")
    A,B = construct_AB_vec(mesh,IOR_dict,k,sparse=sparse)
    N = A.shape[0]

    if A.shape[0]>2000 and not ignore_warning and not sparse:
        raise Exception("A and B matrices are larger than 2000 x 2000 - this may make your system unstable. consider setting sparse=True")
    
    if verbose:
        print("solving...")

    if not sparse:
        _w,_v = eig(A,B,overwrite_a=True,overwrite_b=True)
        inds = _w.argsort()[::-1]
        w = _w[inds][:Nmax]
        v = _v[:,inds][:,:Nmax]
    else:
        if sparse_solve_mode == "transform":
            C = spsolve(A - est_eigval*B,B.todense())
            w,v = eigs(C,Nmax,which='SR')
            w = est_eigval + 1/w
        elif sparse_solve_mode == "pardiso":
            C = pypardiso.spsolve(A - est_eigval*B,B.todense())
            w,v = eigs(C,Nmax,which='SR') 
            w = est_eigval + 1/w
        else:
            w,v = eigsh(A,k=Nmax,M=B,which='SA',sigma=est_eigval)
    if verbose:
        print("solving complete")
    IORs = [ior[1] for ior in IOR_dict.items()]
    nmin,nmax = min(IORs) , max(IORs)
    mode_count = 0
    for _w,_v in zip(w,v.T):
        if _w<0:
            continue
        ne = np.sqrt(_w/k**2)
        if plot:
            if not (nmin <= ne <= nmax) and target_neff is None:
                print("warning: spurious mode! stopping plotting ... ")
            print("effective index: ",get_eff_index(wl,_w))
            plot_vector_mode(mesh,_v)
        if (nmin <= ne <= nmax) or target_neff is not None:
            mode_count+=1
        else:
            break
    return w,v.T,mode_count

#endregion

#region misc

def get_eff_index(wl,w):
    """ get effective index from wavelength wl and eigenvlaue w """
    k = 2*np.pi/wl
    return np.sqrt(w/k**2)

def compute_diff(tri_idx,mesh,_pinv):
    from wavesolve.shape_funcs import compute_NN
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

#endregion
    
#region plotting
    
def plot_eigenvector(mesh,v,show_mesh = False,ax=None,show=True):
    print("deprecated - switch to plot_scalar_mode() or plot_vector_mode()")
    points = mesh.points
    if ax is None:
        fig,ax = plt.subplots(figsize=(5,5))
    
    ax.set_aspect('equal')
    im = ax.tricontourf(points[:,0],points[:,1],v,levels=60)
    
    if show_mesh:
        plot_mesh(mesh,ax=ax)
    if show:
        plt.show()
    return im

def plot_scalar_mode(mesh,v,show_mesh=False,ax=None):
    """ plot a scalar eigenmode 
    ARGS
        mesh: finite element mesh
        v: an array (column vector) corresponding to an eigenmode
        show_mesh: set True to additionally plot the mesh geometry
        ax: optionally put the plot on a specific matplotlib axis
    """
    points = mesh.points
    show=False
    if ax is None:
        show=True
        fig,ax = plt.subplots(figsize=(5,5))
    
    ax.set_aspect('equal')
    im = ax.tricontourf(points[:,0],points[:,1],v,levels=60)
    
    if show_mesh:
        plot_mesh(mesh,ax=ax)
    if show:
        plt.show()
    return im    

def plot_vector_mode(mesh,v,show_mesh=False,ax=None,arrows=True):
    """ plot a scalar eigenmode 
    ARGS
        mesh: finite element mesh
        v: an array (column vector) corresponding to an eigenmode
        show_mesh: set True to additionally plot the mesh geometry
        ax: optionally put the plot on a specific matplotlib axis
        arrows: whether or not to overplot field arrows
    """
    tris = mesh.cells[1].data

    edge_inds = mesh.edge_indices
    show = False
    if ax is None:
        fig,ax = plt.subplots(1,1)
        ax.set_aspect('equal')
        show = True

    amps = []
    vecs = []
    xps = []
    yps = []

    for tri,edge in zip(tris,edge_inds):
        tripoints = mesh.points[tri]
        centroid = np.mean(tripoints,axis=0)

        vec = LNe0(centroid,tripoints,tri)*v[edge[0]] + LNe1(centroid,tripoints,tri)*v[edge[1]] +LNe2(centroid,tripoints,tri)*v[edge[2]]
        amps.append(np.linalg.norm(vec))
        vecs.append(vec)
        xps.append(centroid[0])
        yps.append(centroid[1])

    vecs = np.array(vecs)

    ax.tricontourf(xps,yps,amps,levels=60)
    if arrows:
        ax.quiver(xps,yps,vecs[:,0],vecs[:,1],color='white')
    if show_mesh:
        plot_mesh(mesh,ax=ax)
    if show:
        plt.show()

#endregion

#region field evaluation

def det(u,v):
    return u[0]*v[1] - u[1]*v[0]

def isinside(v, tri, include_vertex = True):
    ''' checks if the given point is inside the triangle '''
    v0 = tri[0]
    v1 = tri[1] - tri[0]
    v2 = tri[2] - tri[0]

    a = (det(v,v2) - det(v0,v2)) / det(v1,v2)
    b = -(det(v,v1) - det(v0,v1)) / det(v1,v2)

    if include_vertex:
        if (a>=0 and b>=0 and a+b<=1): 
            return True
        else: return False    
    else:
        if (a>0 and b>0 and a+b<1): 
            return True
        else: return False    

def find_triangle(gridpoint, mesh):
    ''' 
    Finds which triangle the point lies in

    Args: 
    gridpoint: [x,y] cartesian coordinates
    mesh: mesh
    
    Output:
    triangle_index: the index of the triangle that the [x,y] point lies in.
                    returns -99 if the point doesn't lie in any triangle.
    '''
    points = mesh.points
    tris = mesh.cells[1].data 
    for i in range(len(tris)):
        tri_points = points[tris[i]]
        if isinside(gridpoint, tri_points[:,:2]):
            return i
    return None

def interpolate_field(gridpoint, index, v, mesh):
    '''
    Finds the field at [x,y] by interpolating the solution found on the triangular mesh

    Args:
    gridpoint: [x,y] cartesian coordinates
    index: the index of the triangle that the point lies in (found from find_triangle function)
    v: the field (solution) to interpolate
    mesh: mesh

    Output:
    interpolated: the interpolated field on the [x, y] point
    '''

    if index*0 != 0: return np.nan
    points = mesh.points
    tris = mesh.cells[1].data
    field_points = v[tris[int(index)]]

    vertices = points[tris[int(index)]][:,:2]
    uvcoord = affine_transform(vertices)(gridpoint)
    interpolated = 0
    for ii in range(6):
        interpolated += get_basis_funcs_affine()[ii](uvcoord[0], uvcoord[1]) * field_points[ii]

    return interpolated 

def get_tri_idxs(mesh,xa,ya):
    tri_idxs = np.zeros((len(xa), len(ya)),dtype=int)
    for i in range(len(xa)):
        for j in range(len(ya)):
            idx = find_triangle([xa[i],ya[j]], mesh) 
            tri_idxs[i][j] = idx if idx is not None else -1
    return tri_idxs

def get_interp_weights(mesh,xa,ya,tri_idxs,order=2):
    assert order in [1,2], "order must be 1 or 2, corresponding to mesh element order"

    if order==2:
        N = 6
    else:
        N = 3

    weights = np.zeros((len(xa),len(ya),N))

    points = mesh.points
    tris = mesh.cells[1].data

    for i in range(len(xa)):
        for j in range(len(ya)):
            for k in range(N):
                if tri_idxs[i,j]==-1:
                    weights[i,j,k] = np.nan
                    continue
                gridpoint = [xa[i],ya[j]]
                vertices = points[tris[tri_idxs[i,j]]][:,:2]
                gridpoint_uv = affine_transform(vertices)(gridpoint)
                if order==2:
                    weights[i,j,k] = get_basis_funcs_affine()[k](gridpoint_uv[0], gridpoint_uv[1])
                else:
                    weights[i,j,k] = get_linear_basis_funcs_affine()[k](gridpoint_uv[0], gridpoint_uv[1])
    return weights

def get_interp_weights_vec(mesh,xa,ya,tri_idxs):
    N = 3
    weights = np.zeros((len(xa),len(ya),N,2)) # weights are vectorial; list dim stores (x,y) component

    points = mesh.points
    tris = mesh.cells[1].data

    for i in range(len(xa)):
        for j in range(len(ya)):
            for k in range(N):
                tri = tri_idxs[i,j]
                if tri==-1:
                    weights[i,j,k] = np.nan
                    continue
                
                gridpoint = [xa[i],ya[j]]
                vertices = points[tris[tri]][:,:2]

                weights[i,j,k,:] = get_edge_linear_basis_funcs_affine()[k](gridpoint,vertices,mesh.cells[1].data[tri]) #* mesh.edge_flips[tri,k]
    return weights

def interpolate(v,mesh,xa,ya,tri_idxs = None,interp_weights = None,meshtree=None,order=2,maxr=0):
    """ interpolates SCALAR eigenmode v, computed on mesh, onto rectangular grid defined by 1D arrays xa and ya.
    ARGS:
        v: eigenvector to interpolate 
        mesh: mesh object corresponding to waveguide geometry
        xa: 1D array of x points for output grid
        ya: 1D array of y points for output grid
        tri_idxs: an array of indices. the first index corresponds to the first triangle containing the first mesh point, etc.
        interp_weights: interpolation weights. these are multiplied against v and summed to get the interpolated field
        mesh_tree: a KDtree representing the mesh triangles. if None, one will be made with construct_meshtree(mesh)
        order: the order of the finite element mesh
        maxr: points more than this distance from the origin are ignored. if 0, all points are assumed to lie inside the mesh
    RETURNS:
        the mode v interpolated over the rectangular grid (xa,ya)
    """

    if meshtree is None:
        meshtree = construct_meshtree(mesh)

    tri_idxs = get_tri_idxs_KDtree(mesh,meshtree,xa,ya,maxr=maxr) if tri_idxs is None else tri_idxs
    interp_weights = get_interp_weights(mesh,xa,ya,tri_idxs,order=order) if interp_weights is None else interp_weights

    tris = mesh.cells[1].data
    field_points = v[tris[tri_idxs]]

    return np.sum(field_points*interp_weights,axis=2)

def interpolate_vec(v,mesh,xa,ya,tri_idxs = None,interp_weights = None,meshtree=None,maxr=0):
    """ interpolates the VECTOR eigenmode v, computed on mesh, onto rectangular grid defined by 1D arrays xa and ya.
    ARGS:
        v: eigenvector to interpolate 
        mesh: mesh object corresponding to waveguide geometry
        xa: 1D array of x points for output grid
        ya: 1D array of y points for output grid
        tri_idxs: an array of indices. the first index corresponds to the first triangle containing the first mesh point, etc.
        interp_weights: interpolation weights. these are multiplied against v and summed to get the interpolated field
        mesh_tree: a KDtree representing the mesh triangles. if None, one will be made with construct_meshtree(mesh)
        order: the order of the finite element mesh
        maxr: points more than this distance from the origin are ignored. if 0, all points are assumed to lie inside the mesh
    RETURNS:
        the mode v interpolated over the rectangular grid (xa,ya)
    """
    edge_indices = mesh.edge_indices
    edges = mesh.cells[0].data # for this to work, need to update mesh with get_unique_edges()
    Ntt = len(edges)
    vtt = v[:Ntt]

    if meshtree is None:
        meshtree = construct_meshtree(mesh)

    tri_idxs = get_tri_idxs_KDtree(mesh,meshtree,xa,ya,maxr=maxr) if tri_idxs is None else tri_idxs
    interp_weights = get_interp_weights_vec(mesh,xa,ya,tri_idxs) if interp_weights is None else interp_weights

    field_points = vtt[edge_indices[tri_idxs]]

    return np.sum(field_points[:,:,:,None]*interp_weights,axis=2)

def unstructured_interpolate(v,inmesh,point,inmeshtree=None,max_tries = 10):
    """ interpolate the field v evaluated on the points of inmesh to compute the value at an arbitrary [x,y] point. 
        the function is sped up drastically through the use of a KDtree, which must be supplied for the inmeshtree arg.
        to get a KDtree for a mesh, run mesher.construct_meshtree(mesh).

        this function uses a KDtree to progressively search through triangles in the mesh, ordering the search in terms
        of how close the interpolation point is to each triangle centroid. after searching through max_tries triangles,
        the algorithm will assume the point is outside the mesh and return 0.
    
    ARGS:
        v: the mode field you want to resample
        inmesh: the mesh that v was computed on
        point: the [x,y] point you want to interpolate to
        inmeshtree: the KDtree for inmesh. if none, it is auto-computed
        max_tries: number of searches in KDtree done before we assume that the point is outside the mesh
    RETURNS: 
        the interpolated value of v at point. returns nan if outside the mesh (as per max_tries criterion)
    """
    tryno = 0
    if inmeshtree is None:
        inmeshtree = construct_meshtree(inmesh)

    while tryno < max_tries:
        tri_idx = inmeshtree.query(point,[tryno+1])[1][0]
        tri_idxs = inmesh.cells[1].data[tri_idx]
        tripoints = inmesh.points[tri_idxs,:2]
        if not isinside(point,tripoints):
            tryno += 1 
            continue
        uv = apply_affine_transform(tripoints,point)
        Ns = evaluate_basis_funcs(*uv)
        return np.sum(Ns*v[tri_idxs])
    return np.nan

def mesh_interpolate(v,inmesh,outmesh,inmeshtree=None,max_tries = 10):
    """ interpolate the field v, sampled on inmesh.points, to outmesh.points 
    ARGS:
        v: the mode field you want to resample
        inmesh: the mesh that v was computed on
        outmesh: the mesh that v will be resampled on
        inmeshtree: the KDtree for inmesh. if none, it is auto-computed
        max_tries: number of searches in KDtree done before we assume that the point is outside the mesh
    RETURNS:
        the resampled field 
    """
    
    if inmeshtree is None:
        inmeshtree = construct_meshtree(inmesh)

    vout = np.zeros(outmesh.points.shape[0])
    for i,point in enumerate(outmesh.points):
        vout[i] = unstructured_interpolate(v,inmesh,point[:2],inmeshtree,max_tries)
    return vout

def get_mesh_interpolate_matrices(inmesh,outmesh,inmeshtree,boundary_val=0):
    """get the index and weight matrices that interpolate between inmesh and outmesh"""
    weight_matrix = np.empty((outmesh.points.shape[0],6))
    idx_matrix = np.zeros((outmesh.points.shape[0]),dtype=int)
    for i,point in enumerate(outmesh.points):
        idx = find_triangle_KDtree(point,inmesh,inmeshtree,max_tries=-1)
        tridx = inmesh.cells[1].data[idx]
        if idx is not None:
            tripts = inmesh.points[tridx,:2]
            uv = affine_transform(tripts)(point)
            Ns = evaluate_basis_funcs(*uv)
            idx_matrix[i] = tridx
            weight_matrix[i,:] = Ns
        else:
            weight_matrix[i,:] = boundary_val
    
    return idx_matrix,weight_matrix

def find_triangle_KDtree(point,mesh,meshtree,max_tries=-1,maxr = 0):
    if maxr > 0:
        if np.sqrt(point[0]**2+point[1]**2)>maxr:
            return -1
    tryno = 0
    if max_tries == -1:
        max_tries = mesh.points.shape[0]-1
    while tryno < max_tries:
        tri_idx = meshtree.query(point,[tryno+1])[1][0] # i know this is inefficient ...
        tri_idxs = mesh.cells[1].data[tri_idx]
        tripoints = mesh.points[tri_idxs,:2]
        if not isinside(point,tripoints):
            tryno += 1 
            continue
        return tri_idx
    return None

def get_tri_idxs_KDtree(mesh,meshtree,xa,ya,maxr=0):
    tri_idxs = np.zeros((len(xa), len(ya)),dtype=int)
    for i in range(len(xa)):
        for j in range(len(ya)):
            idx = find_triangle_KDtree([xa[i],ya[j]],mesh,meshtree,maxr=maxr) 
            tri_idxs[i][j] = idx if idx is not None else -1
    return tri_idxs

#endregion