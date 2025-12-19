----------------------------------------
``wavesolve``: a waveguide mode solver 
----------------------------------------

``wavesolve`` is a lightweight package which can solve for the scalar and vector eigenmodes of waveguides using the finite-element method (specifically, Galerkin's method - read more `here <https://github.com/jw-lin/wavesolve/blob/julia/finite_element_method_notes.pdf/>`_). This particular branch offloads some calculations to Julia while maintaining a completely Pythonic interface, and at some point will replace the original pure Python branch.

The general workflow for ``wavesolve`` is the following:

1. Define a waveguide cross-section, by setting refractive index values and material boundaries
2. Generate a finite element mesh for the waveguide
3. Solve the generalized eigenproblem on the mesh, giving the eigenmodes of the waveguide.
4. Plot and resample the eigenmodes as needed.

To implement the above, ``wavesolve`` provides the following submodules:

1. ``waveguide``: a set of classes which define waveguides and generate meshes.
2. ``FEsolver``: a set of functions for eigenmode solving, plotting, and resampling.

.. note::

   ``wavesolve`` is a work in progress. There may be bugs.

^^^^^^^^^^^^^^^^^
contents
^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   installation
   quickstart
   examples
   advanced
   reference
