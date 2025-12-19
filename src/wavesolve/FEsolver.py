
import numpy as np
import os,wavesolve
from wavesolve.waveguide import plot_mesh
import matplotlib.pyplot as plt

from juliacall import Main as jl
from juliacall import Pkg as jlPkg
jlPkg.activate(os.path.dirname(wavesolve.__file__)+"/FEsolver")
jl.seval("using FEsolver")

# region solving
def count_modes(w,wl,IOR_dict):
    """ count the number of guided (non-spurious) modes. This function only works well for weakly guiding waveguides!
        
    ARGS:
        w (array): array of eigenvalues corresponding to the modes, assumed decreasing.
        wl (float): wavelength
        IOR_dict (dict): a dictionary assigning different named regions of the mesh different refractive index values
    """
    nmin,nmax = min(IOR_dict.values()),max(IOR_dict.values()) 
    mode_count = 0
    for _w in w:
        if _w<0:
            continue
        ne = get_eff_index(wl,_w)
        if (nmin <= ne <= nmax):
            mode_count+=1
        else:
            break

    return mode_count

def solve_waveguide(mesh,wl,IOR_dict,Nmax=6,target_neff=None,solve_mode="sparse"):
    """ given a mesh, propagation wavelength, and refractive index dictionary, solve for the scalar modes. 
    
    ARGS: 
        mesh: mesh object corresponding to waveguide geometry
        wl (float): wavelength, defined in the same units as mesh point positions
        IOR_dict (dict): a dictionary assigning different named regions of the mesh different refractive index values
        Nmax (int): return only the <Nmax> largest eigenvalue/eigenvector pairs
        target_neff (float or None): search for modes with indices close to but below this value. if None, target_neff is set to the maximum index in the guide.
        solve_mode (str): use 'sparse' or 'dense' matrices when solving
    
    RETURNS:
        (tuple) : a tuple containing:
            
            - eigvals: array of eigenvalues, descending order
            - eigvecs: array of corresponding eigenvectors (waveguide modes)
    """

    assert solve_mode in ["sparse","dense"], "solve_mode must be 'sparse' or 'dense' (default 'sparse')."

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
    w,v = jl.FEsolver.solve_waveguide(points,tris,IORs,k**2,est_eigval,Nmax,order=order,solve_mode=solve_mode)
    w,v = np.array(w),np.array(v).T

    return w,v

def solve_waveguide_vec(mesh,wl,IOR_dict,Nmax=6,target_neff=None,solve_mode="mixed"):
    """ given a mesh, propagation wavelength, and refractive index dictionary, solve for the vector modes.
    
    ARGS: 
        mesh: mesh object corresponding to waveguide geometry
        wl (float): wavelength, defined in the same units as mesh point positions
        IOR_dict (dict): a dictionary assigning different named regions of the mesh different refractive index values
        Nmax (int): return only the <Nmax> largest eigenvalue/eigenvector pairs
        target_neff (float or None): search for modes with indices close to but below this value. if None, target_neff is set to the maximum index in the guide.
        solve_mode (str): 'mixed' combines a sparse linear solve and a dense eigensolve. 'dense' uses dense matrices for everything; default 'mixed'

    RETURNS:
        (tuple) : a tuple containing:
            
            - eigvals: array of eigenvalues, descending order
            - eigvecs: array of corresponding eigenvectors (waveguide modes)
    """
    assert get_mesh_order(mesh) == 1, "only order 1 meshes supported for vectorial solving"
    assert solve_mode in ["mixed","dense"], "solve_mode must be 'mixed' or 'dense' (default 'mixed')."

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

    w,v = jl.FEsolver.solve_waveguide_vec(points,tris,edges,IORs,k**2,Nedges,est_eigval,Nmax,solve_mode=solve_mode)
    w,v = np.array(w),np.array(v).T

    return w,v

#endregion

#region plotting

def plot_vector_field(mesh,field,show_mesh=False,ax=None,arrows=True,res=100,bounds=None):
    """ plot a vector field - transverse and longitudinal components.

    ARGS:
        mesh: finite element mesh
        field (array): a 1D array corresponding to a vector-valued field (e.g. eigenmode)
                which encodes both transverse and longitudinal components
        show_mesh (bool): set True to additionally plot the mesh geometry
        ax (matplotlib.axes or None): optionally put the plots on mpl axis. if one axis is given, only the transverse part is plotted.
                if (ax0,ax1) is given, the longitudinal component is also plotted. if None, axes are made for both.
        arrows (bool): whether or not to overplot field arrows
        res (int): grid resolution, make tuple to set xres and yres separately
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
        im1 = ax[1].imshow(np.real(vi_z).T,extent=(xa[0],xa[-1],ya[0],ya[-1]),origin="lower")
        plt.colorbar(im1,ax=ax[1],fraction=0.046, pad=0.04)
        if show_mesh:
            plot_mesh(mesh,ax=ax[1])
            ax[1].set_xlim(bounds[0],bounds[1])
            ax[1].set_ylim(bounds[2],bounds[3])
    if show:
        plt.show()


def plot_scalar_field(mesh,field,show_mesh=False,ax=None,res=100,bounds=None):
    """ plot a scalar field

    ARGS:
        mesh: finite element mesh
        field (array): a 1D array corresponding to a scalar-valued field (e.g. eigenmode)
        show_mesh (bool): set True to additionally plot the mesh geometry
        ax (matplotlib.axes or None): optionally put the plots on an axis. if one axis is given, only the transverse part is plotted.
            if (ax0,ax1) is given, the longitudinal component is also plotted. if None, axes are made for both.
        res (int): grid resolution, make tuple to set xres and yres separately
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

#endregion

# region interpolation

def isvectorial(field,mesh):
    """ determine if field is a vector or scalar field 
    
    ARGS:
        field (array): electric field evaluated on a finite element mesh
        mesh: the finite element mesh
    """
    if len(field) == mesh.points.shape[0]:
        return False
    return True

def split_vector_field(field,mesh):
    """split a vectorial field into its transverse and longitudinal parts."""
    assert isvectorial(field,mesh), "field is not vectorial"
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

def evaluate(point,field,mesh):
    """ evaluate a field sampled over a finite element mesh at a given point.
    
    ARGS:
        point: an [x,y] point, or an Nx2 array of points
        field: a real-valued field represented on a finite-element mesh
        mesh: the FE mesh, which has been sorted with sort_mesh().

    RETURNS:
        (float or array): the field evaluated at point(s)
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
    # sort the mesh if necessary
    if not hasattr(newmesh,"tree"):
        sort_mesh(newmesh)
    return evaluate(newmesh.points,field,mesh)

def get_mesh_bounds(mesh):
    return np.array(mesh.tree._idxtree.bbox)

#endregion

#region misc

def sort_mesh(mesh):
    """ create a bounding volume hierarchy tree for ``mesh`` and pass it into to ``mesh.tree``.
    this is required for evaluations of fields on meshes. you usually don't need to 
    call this manually.
    
    ARGS:
        mesh: finite element mesh
    """

    mesh.tree = create_tree_from_mesh(mesh)
    return mesh

def get_mesh_order(mesh):
    return int(mesh.cells[1].data.shape[1]/3)

def dec_2D(v,decim_factor=4):
    return v[::decim_factor,::decim_factor]

def get_eff_index(wl,w):
    """ get effective refractive index from wavelength wl and eigenvlaue w 
    
    ARGS:
        wl (float): wavelength
        w (float or array): eigenvalue(s)
    """
    k = 2*np.pi/wl
    return np.sqrt(w/k**2)

#endregion