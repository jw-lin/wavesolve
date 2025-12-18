Advanced usage
==============

.. contents::
    :local:
    :depth: 1

Making a custom ``Waveguide``
-----------------------------

In this section we discuss how to build a custom ``Waveguide`` from the bottom-up from basic geometric shapes.

``Prim2D`` class
^^^^^^^^^^^^^^^^

The basic parent class that represents a refractive index geometry is a ``Prim2D`` (2D primitive). Each ``Prim2D`` stores a refractive index value and an array of :math:`(x,y)` points bounding the physical region which contains that index value. A ``Prim2D`` is initialized as ::

    prim2D = Prim2D(n,label,points)

**Arguments**

1.  ``n`` : the refractive index of the primitive. 

2. ``label`` : a string identifier to attach to the primitive (e.g. "core" or "cladding").

3. ``points`` (optional): a 2-element array of :math:`(x,y)` points corresponding to a material interface. The first and last point are automatically connected. Default is ``None``, since ``points`` can also be generated with functional dependence, as mentioned next.  

To make custom geometries, define a subclass that inherits from ``Prim2D`` and generate ``points`` according to specific rules. These subclasses should implement their own ``make_points()`` functions, which at minimum should take in some set of arguments (like radius, for a circle primitive) and return the corresponding point array **in counterclockwise order**. Subclasses should also implement the following function:

``boundary_dist(x,y)`` : compute the signed minimum distance between the point :math:`(x,y)` and the primitive boundary, which is negative value if the point is inside the boundary. 

See the ``Circle``, ``Ellipse``, and ``Rectangle`` classes for examples. 

More complicated primitives can be created through ``waveguide.Prim2DUnion``, representing a list of intersecting primitives, and ``waveguide.Prim2DArray``, representing a list of non-intersecting primitives.

.. note::
    ``points`` should always be ordered counterclockwise!

``Waveguide`` construction
^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``Waveguide`` stores a (potentially nested) list of ``Prim2D`` objects, which I will call ``prim2Dgroups``. The refractive index profile of each top-level element in ``prim2Dgroups`` is overwritten by the immediate next element. An element in ``prim2Dgroups`` can also be a list of ``Prim2D`` objects; in this case all elements in the sublist override the refractive index profile of the previous element. [1]_ So, an optical fiber might store its primitives as ``[cladding,core]``, where ``cladding`` and ``core`` are ``Circle(Prim2D)`` objects corresponding to the cladding and core regions; a multicore fiber could have the structure ``[cladding,[core1,core2,...]]``. A ``Waveguide`` is initialized as ::

    wvg = Waveguide(prim2Dgroups)

**Arguments** 

1.  ``prim2Dgroups`` : the potentially nested list of 2D primitives mentioned above.

All ``Waveguide`` objects take multiple ``Prim2D`` objects and arrange them as desired in this list. The ``Waveguide`` class is mainly used to generate meshes, which can be tuned via ``Waveguide`` class attributes, as will be shown in the next section. For example ``Waveguide`` subclasses, check out the ``CircularFiber``, ``EllipticalFiber``, and various photonic-lantern-related classes.

.. [1] ``Prim2D`` objects in the same sublist should never intersect. However, they can have different ``label``\ s and refractive indices. If you want to model intersecting primitives, you should make a ``Prim2DUnion``. See :doc:`PL3`  for an example.

Once the waveguide is initialized, you can get the refractive index dictionary of the waveguide with ``Waveguide.assign_IOR()``. This is what passes material information to the eigensolver.

A useful plotting method for checking your work is ``Waveguide.plot_boundaries()``, which plots all material boundaries in waveguide.

Mesh generation
---------------

Mesh generation is handled by the function ::

    Waveguide.make_mesh(algo=6,order=2,adaptive=True)

**Optional Arguments**

1. ``algo`` : the mesh generation algorithm used by Gmsh. The default ``6`` is good starting point.
2. ``order`` : the order of the finite elements. Order ``1`` (linear triangle) and ``2`` (quadratic triangle) are implemented.
3. ``adaptive``: whether or not to refine mesh element sizes so that the mesh is most refined only at primitive boundaries.

Adaptive meshing is slower but tends to give more accurate meshes and ultimately faster solve times. In adaptive meshing mode, the target triangle size at a given :math:`(x,y)` point is computed as 

.. math::
    
    {\rm target \, size} = d_0\left(1+ \dfrac{s \, d(x,y)}{d_0} \right)^p

where :math:`d_0` is the "default" mesh size corresponding to the resolution of the ``Prim2D`` boundary, :math:`d(x,y)` is the distance between the point :math:`(x,y)` and the primitive's boundary, and :math:`s` and :math:`p` are variables; higher values mean that mesh size will increase more rapidly away from the boundary. For multiple primitives, a target size is computed for each and the minimum size is taken. Then the target size is clipped between a minimum and maximum allowed value. The parameter values for the adaptive scheme are set through the following ``Waveguide`` class attributes:

* ``mesh_dist_scale`` : :math:`s`, the mesh boundary refinement linear distance scaling. Default ``0.25``
* ``mesh_dist_power`` : :math:`p`, mesh boundary refinement power scaling. Default ``1`` (linear)
* ``min_mesh_size`` : minimum allowed mesh size, default ``0.1``
* ``max_mesh_size`` : maximum allowed mesh size, default ``10``

Users can also specify a target mesh size, and toggle boundary refinement on a per-primitive basis. This is done through the following ``Prim2D`` attributes: 

* ``mesh_size`` : target mesh size within the boundary of the primitive (otherwise the mesh size is set by the scheme described above.)
* ``skip_refinement`` : whether or not mesh refinement should be applied at the primitive boundary. The outer boundary of the entire mesh should have this set to ``True``; default ``False``.

To view meshes, the ``Waveguide`` class implements ::

    Waveguide.plot_mesh(mesh=None,IOR_dict=None)

**Optional arguments**

1. ``mesh``: the mesh to plot. If ``None``, one is generated using default values through ``make_mesh()``.
2. ``IOR_dict``: dictionary of refractive index values. If ``None``, one is generated through ``assign_IOR()``.

Scalar modes
------------

Solve
^^^^^

Scalar eigenmodes can be solved on **both** order 1 and order 2 meshes. In either case, we use ::

    eigvals, eigvecs = solve_waveguide(mesh,wl,IOR_dict,Nmax=6,target_neff=None)

**Arguments**

1. ``mesh``: mesh object corresponding to waveguide geometry
2. ``wl``: wavelength, defined in the same units as mesh point positions
3. ``IOR_dict``: a dictionary assigning different named regions of the mesh different refractive index values
4. ``Nmax``: (optional) return only the `Nmax` largest eigenvalue/eigenvector pairs; default 6.
5. ``target_neff``: (optional) upper bound on expected effective index of the mode you want to solve for; autocomputed if `None`.

Note that the autocomputed ``target_neff`` will typically be too high for waveguides which guide light in air or other low-index materials, e.g. :doc:`hollowPCF`.

To count the number of guided modes, use ::

    num_modes = FEsolver.count_modes(w,wl,IOR_dict)

This function will only work for waveguides with weak index contrast.

Plot
^^^^

To plot a scalar field, e.g. a scalar eigenmode (which is a row of the ``eigvecs`` matrix above), use ::

    plot_scalar_field(mesh,field,show_mesh=False,ax=None,res=100,bounds=None)

**Arguments**

1.  ``mesh``: finite element mesh
2.  ``field``: an array (column vector) corresponding to an eigenmode
3.  ``show_mesh``: (optional) set True to additionally plot the mesh geometry
4.  ``ax``: (optional) a matplotlib axis to plot on
5.  ``res``: (optional) grid resolution, make tuple to set xres and yres separately
6.  ``bounds``: (optional) 4-element array ``[xmin,xmax,ymin,ymax]``, setting plot bounds. If ``None``, use the mesh bounds.

Evaluate
^^^^^^^^

To evaluate a finite-element field at an arbitrary point, ``wavesolve.FEsolver`` provides two functions. First, there is ::

    field_amp = evaluate(point,field,mesh)

**Arguments**

1. ``point`` : a 2-element array ``[x,y]`` *or* an :math:`N\times 2` array of such points.
2. ``field`` : the electric field (e.g. an eigenmode) to be evaluated
3. ``mesh`` : the finite element mesh on which ``field`` is defined.

**Returns**

1. ``field_amp`` : the field amplitude at ``point``, or an array of field amplitudes for each point in ``point``.

To evaluate a ``field`` on a rectangular grid, use ::

    field_amp_grid = evaluate(xpoints,ypoints,field,mesh)

1. ``xpoints`` : a 1D array of :math:`x` coordinates to evaluate on
2. ``ypoints`` : a 1D array of :math:`y` coordinates to evaluate on
3. ``field`` : the electric field (e.g. an eigenmode) to be evaluated
4. ``mesh`` : the finite element mesh on which ``field`` is defined.

**Returns**

1. ``field_amp_grid`` : a 2D array of field amplitudes corresponding to the grid defined by ``pointsx`` and ``pointsy``.

.. note::
    Grid evaluation functions (and ``wavesolve`` in general) treat the **first** axis of a matrix as :math:`x` and the **second** as :math:`y`. Therefore, if you use a plotting function such as ``matplotlib.pyplot.imshow``, you need to transpose the field data (and set ``origin="lower"``).

Vector modes
------------

``wavesolve`` divides field vectors into transverse and longitudinal components. The longitudinal component is treated in the same way as a scalar mode, which has a defined value at every node position. In contrast, the transverse electric field is specified by assigning a value for the tangential electric field to every mesh edge. If :math:`N_e` is the number of mesh edges and :math:`N_n` is the number of nodes, a vector mode is an :math:`N_e+N_n` array whose first :math:`N_e` elements determine the transverse field and whose last :math:`N_n` components determine the longitudinal field.

.. note:: ``wavesolve`` performs a change of variables to treat both the transverse and longitudinal field components as purely real, when physically they should be 90 degrees out of phase (i.e. the longitudinal field should be purely imaginary). See `*Lee et al. 1991* <https://ieeexplore.ieee.org/document/85399>`_ for more details.


Solve
^^^^^

The syntax for modesolving is basically the same as the in the scalar case, except that all input meshes must be order 1. Use ::

    eigvals, eigvecs = solve_waveguide_vec(mesh,wl,IOR_dict,Nmax=6,target_neff=None)

**Arguments**

1. ``mesh``: an order 1 mesh object corresponding to waveguide geometry
2. ``wl``: wavelength, defined in the same units as mesh point positions
3. ``IOR_dict``: a dictionary assigning different named regions of the mesh different refractive index values
4. ``Nmax``: (optional) return only the ``Nmax`` largest eigenvalue/eigenvector pairs
5. ``target_neff``: (optional) expected effective index of the mode you want to solve for; autocomputed if ``None``.

**Returns**

1. ``eigvals``: array of eigenvalues, descending order
2. ``eigvecs``: array of vector-valued eigenmodes

The eigenvectors and eigenvalues returned by ``wavesolve`` should always be real, but the data is returned as complex just in case.

The default solving mode is to transform the generalized eigenproblem :math:`A x =\lambda B x` is transformed to :math:`Cx = \lambda x` where :math:`C` solves :math:`BC = A`; this system is solved with ``Pardiso.jl``.

To count the guided modes, I generally recommend to distinguish between the physical and nonphysical modes by eye. However, for weak index constrast, ``count_modes()`` will still work.

Plot
^^^^

To plot a vector eigenmode, use ``plot_vector_field()``, which has the same signature as ``plot_scalar_field()`` but produces a quiver plot corresponding to the transverse component of the electric field, and optionally a scalar plot of the longitudinal field as well. ::

    plot_vector_field(mesh,field,show_mesh=False,ax=None,arrows=True,res=100,bounds=None)

**Arguments**

1. ``mesh``: finite element mesh
2. ``field``: a 1D array (column vector) corresponding to a vector-valued field (e.g. eigenmode) which encodes both transverse and longitudinal components
3. ``show_mesh``: (optional) plot the mesh geometry
4. ``ax``: (optional) put the plot(s) on a matplotlib axis. If one axis is given, only the transverse part is plotted. if (ax0,ax1) is given, the longitudinal component is also plotted. If None, axes are made for both.
5. ``arrows`` : (optional) whether or not to overplot field arrows
6. ``res``: (optional) grid resolution, can be a tuple
7. ``bounds``: 4-element array ``[xmin,xmax,ymin,ymax]``, setting plot boundary.

Evaluate
^^^^^^^^

You can use the same functions as in the scalar case. However, due to the vectorial nature of the fields, evaluation functions will return two quantities: the vectorial transverse field, and the scalar longitudinal field. For the transverse field, the last axis of the matrix will have length 2, corresponding to the :math:`x` and :math:`y` components. For example, if ``xpoints`` and ``ypoints`` are :math:`N`-long coordinate arrays, then for ::

    field_amp_xy, field_amp_z = evaluate_grid(xpoints,ypoints,field,mesh)

``field_amp_xy`` will be an :math:`N\times N\times 2` array and ``field_amp_z`` will be an :math:`N\times N` array.