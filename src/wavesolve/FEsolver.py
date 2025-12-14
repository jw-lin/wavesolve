
import numpy as np
import os,wavesolve
from wavesolve.waveguide import plot_mesh
import matplotlib.pyplot as plt

from juliacall import Main as jl
from juliacall import Pkg as jlPkg
jlPkg.activate(os.path.dirname(wavesolve.__file__)+"/FEsolver")
jl.seval("using FEsolver")

# region solving
def construct_AB(mesh,IOR_dict,k):
    points = mesh.points.T[:2,:]
    materials = mesh.cell_sets.keys()
    all_tris = mesh.cells[1].data
    tris = all_tris.T+1
    IORs = np.empty(tris.shape[1],dtype=np.float64)

    counter = 0
    for material in materials:
        nummat = len(mesh.cell_sets[material][1])
        IORs[counter:counter+nummat] = IOR_dict[material]**2
        counter += nummat

    A,B = jl.FEsolver.construct_AB_order2_sparse(points,tris,IORs,k**2)
    return A,B

def construct_AB_vec(mesh,IOR_dict,k):
    points = mesh.points.T[:2,:]
    materials = mesh.cell_sets.keys()
    all_tris = mesh.cells[1].data
    tris = all_tris.T+1
    edges = mesh.edge_indices.T+1
    Nedges = mesh.cells[0].data.shape[0]
    IORs = np.empty(tris.shape[1],dtype=np.float64)
    counter = 0
    for material in materials:
        nummat = len(mesh.cell_sets[material][1])
        IORs[counter:counter+nummat] = IOR_dict[material]**2
        counter += nummat
    A,B = jl.FEsolver.construct_AB_vec(points,tris,edges,IORs,k**2,Nedges)
    return A,B

def solve(A,B,est_eigval,Nmax):
    w,v = jl.FEsolver.solve(A,B,est_eigval,Nmax)
    return np.array(w),np.array(v).T

def solve_waveguide(mesh,wl,IOR_dict,plot=False,Nmax=6,target_neff=None):
    """ given a mesh, propagation wavelength, and refractive index dictionary, solve for the SCALAR modes. 
        Uses order 2 triangular finite elements.
    
    ARGS: 
        mesh: mesh object corresponding to waveguide geometry
        wl: wavelength, defined in the same units as mesh point positions
        IOR_dict: a dictionary assigning different named regions of the mesh different refractive index values
        plot: set True to view eigenmodes
        Nmax: return only the <Nmax> largest eigenvalue/eigenvector pairs
        target_neff: search for modes with indices close to but below this value. if None, target_neff is set to the maximum index in the guide.
    RETURNS:
        w: array of eigenvalues, descending order
        v: array of corresponding eigenvectors (waveguide modes)
        N: number non-spurious (i.e. propagating) waveguide modes
    """

    # sort the mesh if necessary
    if not hasattr(mesh,"tree"):
        sort_mesh(mesh)

    k = 2*np.pi/wl
    points = mesh.tree.points
    materials = mesh.cell_sets.keys()
    tris = mesh.tree.tris

    IORs = np.empty(tris.shape[1],dtype=np.float64)
    counter = 0
    for material in materials:
        nummat = len(mesh.cell_sets[material][1])
        IORs[counter:counter+nummat] = IOR_dict[material]**2
        counter += nummat

    if target_neff is None:
        est_eigval = np.power(k*max(IOR_dict.values()),2)
    else:
        est_eigval = np.power(k*target_neff,2)

    order = get_mesh_order(mesh)
    w,v = jl.FEsolver.solve_waveguide(points,tris,IORs,k**2,est_eigval,Nmax,order=order)
    w,v = np.array(w),np.array(v).T

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
            plot_scalar_field(mesh,_v)
        if (nmin <= ne <= nmax):
            mode_count+=1
        else:
            break
    return w,v,mode_count

def solve_waveguide_vec(mesh,wl,IOR_dict,plot=False,Nmax=6,target_neff=None,solve_mode="transform"):
    """ given a mesh, propagation wavelength, and refractive index dictionary, solve for VECTOR modes, using linear triangles (order 1).
    
    ARGS: 
        mesh: mesh object corresponding to waveguide geometry
        wl: wavelength, defined in the same units as mesh point positions
        IOR_dict: a dictionary assigning different named regions of the mesh different refractive index values
        plot: set True to view eigenmodes
        Nmax: return only the <Nmax> largest eigenvalue/eigenvector pairs
        target_neff: search for modes with indices close to but below this value. if None, target_neff is set to the maximum index in the guide.

    RETURNS:
        w: array of eigenvalues, descending order
        v: array of corresponding eigenvectors (waveguide modes)
        N: number non-spurious (i.e. propagating) waveguide modes
    """
    assert get_mesh_order(mesh) == 1, "only order 1 meshes supported for vectorial solving"

    # sort the mesh if necessary
    if not hasattr(mesh,"tree"):
        sort_mesh(mesh)

    k = 2*np.pi/wl
    points = mesh.tree.points
    materials = mesh.cell_sets.keys()
    tris = mesh.tree.tris
    edges = mesh.tree.edges
    Nedges = mesh.cells[0].data.shape[0]
    IORs = np.empty(tris.shape[1],dtype=np.float64)
    counter = 0
    for material in materials:
        nummat = len(mesh.cell_sets[material][1])
        IORs[counter:counter+nummat] = IOR_dict[material]**2
        counter += nummat

    if target_neff is None:
        est_eigval = np.power(k*max(IOR_dict.values()),2)
    else:
        est_eigval = np.power(k*target_neff,2)

    w,v = jl.FEsolver.solve_waveguide_vec(points,tris,edges,IORs,k**2,Nedges,est_eigval,Nmax,solve_mode)
    w,v = np.array(w),np.array(v).T

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
            plot_vector_field(mesh,_v)
        if (nmin <= ne <= nmax) or target_neff is not None:
            mode_count+=1
        else:
            break

    return w,v,mode_count

#endregion

#region plotting

def dec_2D(v,decim_factor=4):
    return v[::decim_factor,::decim_factor]

def plot_vector_field(mesh,field,show_mesh=False,ax=None,arrows=True,res=100,bounds=None):
    """ plot a vector field - transverse and (modulus of) longitudinal components. 
    ARGS
        mesh: finite element mesh
        field: a 1D array (column vector) corresponding to a vector-valued field (e.g. eigenmode)
            which encodes both transverse and longitudinal components
        show_mesh: set True to additionally plot the mesh geometry
        ax: optionally put the plots on mpl axis. if one axis is given, only the transverse part is plotted.
            if (ax0,ax1) is given, the longitudinal component is also plotted. if None, axes are made for both.
        arrows: whether or not to overplot field arrows
        res: grid resolution, make tuple to set xres and yres separately
        bounds: 4-element array [xmin,xmax,ymin,ymax], setting plot boundary.
    """
    if bounds is None:
        bounds = get_mesh_bounds(mesh)
    if type(res) is tuple:
        xres,yres = res[0],res[1]
    else:
        xres,yres = res,res

    show=False
    if ax is None:
        show=True
        fig,ax = plt.subplots(1,2,sharey=True)
        ax[0].set_aspect('equal')
        ax[1].set_aspect('equal')

    xa = np.linspace(bounds[0],bounds[1],xres)
    ya = np.linspace(bounds[2],bounds[3],yres)
    xg,yg = np.meshgrid(xa,ya)
    vi_t,vi_z = evaluate_grid(xa,ya,field,mesh)

    if not hasattr(ax, '__len__'):
        ax = [ax]

    im0 = ax[0].imshow(np.linalg.norm(vi_t,axis=2).T,extent=bounds,origin="lower")
    if arrows:
        ax[0].quiver(dec_2D(xg),dec_2D(yg),dec_2D(vi_t[:,:,0]),dec_2D(vi_t[:,:,1]),color="white",alpha=0.5,width=5/1000)

    plt.colorbar(im0,ax=ax[0],fraction=0.046, pad=0.04)
    if show_mesh:
        plot_mesh(mesh,ax=ax[0])
        ax[0].set_xlim(bounds[0],bounds[1])
        ax[0].set_ylim(bounds[2],bounds[3])

    if len(ax)>1:
        ax[0].set_title("transverse")
        ax[1].set_title("longitudinal")
        im1 = ax[1].imshow(np.abs(vi_z).T,extent=(xa[0],xa[-1],ya[0],ya[-1]),origin="lower")
        plt.colorbar(im1,ax=ax[1],fraction=0.046, pad=0.04)
        if show_mesh:
            plot_mesh(mesh,ax=ax[1])
            ax[1].set_xlim(bounds[0],bounds[1])
            ax[1].set_ylim(bounds[2],bounds[3])
    if show:
        plt.show()


def plot_scalar_field(mesh,field,show_mesh=False,ax=None,res=100,bounds=None):
    """ plot a scalar field
    ARGS
        mesh: finite element mesh
        field: a 1D array (column vector) corresponding to a scalar-valued field (e.g. eigenmode)
        show_mesh: set True to additionally plot the mesh geometry
        ax: optionally put the plots on mpl axis. if one axis is given, only the transverse part is plotted.
            if (ax0,ax1) is given, the longitudinal component is also plotted. if None, axes are made for both.
        res: grid resolution, make tuple to set xres and yres separately
        bounds: 4-element array [xmin,xmax,ymin,ymax], setting plot boundary.
    """

    if bounds is None:
        bounds = get_mesh_bounds(mesh)
    if type(res) is tuple:
        xres,yres = res[0],res[1]
    else:
        xres,yres = res,res

    show=False
    if ax is None:
        show=True
        fig,ax = plt.subplots(1,1)
    xa = np.linspace(bounds[0],bounds[1],xres)
    ya = np.linspace(bounds[2],bounds[3],yres)
    vi = evaluate_grid(xa,ya,field,mesh)
    im = ax.imshow(vi.T,origin="lower",extent = bounds)
    plt.colorbar(im,ax=ax,fraction=0.046, pad=0.04)
    if show_mesh:
        plot_mesh(mesh,ax=ax)
        ax.set_xlim(bounds[0],bounds[1])
        ax.set_ylim(bounds[2],bounds[3])
    if show:
        plt.show()

def plot_scalar_mode_old(mesh,v,show_mesh=False,ax=None):
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

def plot_vector_mode_old(mesh,v,show_mesh=False,ax=None,arrows=True):
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
        tripoints = mesh.points[tri].T
        centroid = np.mean(tripoints,axis=1)

        vec = jl.FEsolver.LNe0(centroid,tripoints,tri)*v[edge[0]] + jl.FEsolver.LNe1(centroid,tripoints,tri)*v[edge[1]] + jl.FEsolver.LNe2(centroid,tripoints,tri)*v[edge[2]]
        amps.append(np.linalg.norm(vec))
        vecs.append(vec)
        xps.append(centroid[0])
        yps.append(centroid[1])

    vecs = np.array(vecs)

    im = ax.tricontourf(xps,yps,amps,levels=60)
    plt.colorbar(im,ax=ax)
    if arrows:
        ax.quiver(xps,yps,vecs[:,0],vecs[:,1],color='white')
    if show_mesh:
        plot_mesh(mesh,ax=ax)
    if show:
        plt.show()

#endregion

# region interpolation

def isvectorial(field,mesh):
    """ determine if field is a vector or scalar field """
    if len(field) == mesh.points.shape[0]:
        return False
    return True

def split_vector_field(field,mesh):
    """split a vectorial field into its transverse and longitudinal parts."""
    assert isvectorial(field), "field is not vectorial"
    Ne = mesh.num_edges
    return field[:Ne],field[Ne:]

def create_tree(points,tris,edges):
    """ from an array of mesh points and an index array of (quadratic) triangle connections, 
    construct a bounding volume hierarchy (BVH) tree, which will be used to evaluate fields
    define on the mesh nodes.
    
    ARGS:
        points: an array of the (x,y) positions of the mesh nodes, dimension N x 2 for N nodes.
        connections: an array containing each triangle in the mesh; each triangle is represented
                     as 6 indices, corresponding to 6 points
    RETURNS:
        bvhtree: the BVH tree for the given mesh points and connections.
    """
    return jl.FEsolver.construct_tritree(points[:,:2].T,tris.T+1,edges.T+1)

def create_tree_from_mesh(mesh):
    """ create a BVH tree directly from a finite element mesh object.
    
    ARGS:
        mesh: a meshio object representing a finite element mesh
    RETURNS:
        bvhtree: the BVH tree for the given mesh points and connections.
    """
    return jl.FEsolver.construct_tritree(mesh.points[:,:2].T,mesh.cells[1].data.T+1,mesh.edge_indices.T+1)

def sort_mesh(mesh):
    """ create a BVH tree for ``mesh`` and pass it into to ``mesh.tree`` """

    mesh.tree = create_tree_from_mesh(mesh)
    return mesh

def get_mesh_order(mesh):
    return int(mesh.cells[1].data.shape[1]/3)

def query(point,mesh):
    """ find the index of the triangle in the mesh that contains the given point. 
    
    ARGS:
        point: an array [x,y] corresponding to the query point
        mesh: the FE mesh, which is assumed to be sorted with sort_mesh().
    RETURNS:
        (int): the index of the triangle containing point in the mesh.
    """

    jl_idx = jl.FEsolver.query(point,mesh.tree)
    return jl_idx-1

def evaluate(point,field,mesh):
    """ evaluate a field sampled over a finite element mesh at a given point.
    
    ARGS:
        point: an [x,y] point, or an Nx2 array of points
        field: a real-valued field represented on a finite-element mesh
        mesh: the FE mesh, which has been sorted with sort_mesh().

    RETURNS:
        (float or vector): the field evaluated at point(s)
    """
    if isvectorial(field,mesh):
        if point.ndim == 2:
            p = point[:,:2].T
        else:
            p = point
        Ne = mesh.num_edges
        return np.array(jl.FEsolver.evaluate_vec(p,field[:Ne],mesh.tree)).T , np.array(jl.FEsolver.evaluate(p,field[Ne:],mesh.tree,order=1))
    else:
        order = get_mesh_order(mesh)
        if point.ndim == 2:
            return np.array(jl.FEsolver.evaluate(point[:,:2].T,field,mesh.tree,order=order))
        return np.array(jl.FEsolver.evaluate(point,field,mesh.tree,order=order))

def evaluate_grid(pointsx,pointsy,field,mesh):
    """ evaluate a field defined over a finite element mesh on a cartesian grid.
    
    ARGS:
        pointsx: a 1D array of x points
        pointsy: a 1D array of y points
        field: a real-valued field represented on a finite-element mesh
        mesh: the FE mesh, which has been sorted with sort_mesh().

    RETURNS:
        (array): a 2D or 3D array corresponding to field (scalar or vectorial), evaluated on the grid.
    """
    if isvectorial(field,mesh):
        Ne = mesh.num_edges
        return np.array(jl.FEsolver.evaluate_vec(pointsx,pointsy,field[:Ne],mesh.tree)).transpose(1,2,0), np.array(jl.FEsolver.evaluate(pointsx,pointsy,field[Ne:],mesh.tree,order=1))
    else:
        order = get_mesh_order(mesh)
        return np.array(jl.FEsolver.evaluate(pointsx,pointsy,field,mesh.tree,order=order))

def resample(field,mesh,newmesh):
    """ resample a finite element field onto a new mesh
    
    ARGS: 
        field: the finite element field to be resampled.
        mesh: the finite element mesh on which <field> is defined.
        newmesh: the new finite element mesh <field> should be sampled on.
    """
    return evaluate(newmesh.points,field,mesh.tree)

def get_mesh_bounds(mesh):
    return np.array(mesh.tree._idxtree.bbox)

#endregion

#region misc

def get_eff_index(wl,w):
    """ get effective index from wavelength wl and eigenvlaue w """
    k = 2*np.pi/wl
    return np.sqrt(w/k**2)

#endregion