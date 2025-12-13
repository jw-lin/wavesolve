from juliacall import Main as jl
from juliacall import Pkg as jlPkg
import numpy as np
from scipy.sparse import csc_matrix
import os,wavesolve

jlPkg.activate(os.path.dirname(wavesolve.__file__)+"/FEsolver")
jl.seval("using FEsolver")

def construct_AB(mesh,IOR_dict,k):
    points = mesh.points.T[:2,:]
    materials = mesh.cell_sets.keys()
    all_tris = mesh.cells[1].data
    tris = all_tris.T+1
    IORs = np.empty(points.shape[1],dtype=np.float64)

    counter = 0
    for material in materials:
        nummat = len(mesh.cell_sets[material][1])
        IORs[counter:counter+nummat] = IOR_dict[material]**2
        counter += nummat

    A,B = jl.FEsolver.construct_AB_order2_sparse(points,tris,IORs,k**2)
    return A,B

def construct_AB_py(mesh,IOR_dict,k):
    points = mesh.points.T[:2,:]
    materials = mesh.cell_sets.keys()
    all_tris = mesh.cells[1].data
    tris = all_tris.T+1
    IORs = np.empty(points.shape[1],dtype=np.float64)

    counter = 0
    for material in materials:
        nummat = len(mesh.cell_sets[material][1])
        IORs[counter:counter+nummat] = IOR_dict[material]**2
        counter += nummat

    Is,Js,Avals,Bvals = jl.FEsolver.compute_AB_vals(points,tris,IORs,k**2)
    Is = np.array(Is)-1
    Js = np.array(Js)-1
    A = csc_matrix((Avals,(Is,Js)))
    B = csc_matrix((Bvals,(Is,Js)))
    return A,B

def compute_NN_dNdN(tri):
    return jl.FEsolver.compute_NN_dNdN(tri)

def solve(A,B,est_eigval,Nmax):
    w,v = jl.FEsolver.solve(A,B,est_eigval,Nmax)
    return np.array(w),np.array(v).T

def solve_waveguide(mesh,wl,IOR_dict,Nmax=6,target_neff=None):
    k = 2*np.pi/wl
    points = mesh.points.T[:2,:]
    materials = mesh.cell_sets.keys()
    #all_tris = [mesh.cells[1].data[tuple(mesh.cell_sets[material])][0,:,0,:] for material in materials]
    #all_tris_v = np.vstack(all_tris)
    #tris = all_tris_v.T+1
    tris = mesh.cells[1].data.T+1

    IORs = np.empty(points.shape[1],dtype=np.float64)
    counter = 0
    for material in materials:
        nummat = len(mesh.cell_sets[material][1])
        IORs[counter:counter+nummat] = IOR_dict[material]**2
        counter += nummat

    if target_neff is None:
        est_eigval = np.power(k*max(IOR_dict.values()),2)
    else:
        est_eigval = np.power(k*target_neff,2)

    w,v = jl.FEsolver.solve_waveguide(points,tris,IORs,k**2,est_eigval,Nmax)
    w,v = np.array(w),np.array(v).T

    IORs = [ior[1] for ior in IOR_dict.items()]
    nmin,nmax = min(IORs) , max(IORs)
    mode_count = 0
    for _w,_v in zip(w[::-1],v.T[::-1]):
        if _w<0:
            continue
        ne = np.sqrt(_w/k**2)
        if (nmin <= ne <= nmax):
            mode_count+=1
        else:
            break

    return w,v,mode_count