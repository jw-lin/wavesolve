wavesolve.waveguide
-------------------
.. contents::
    :local:
    :depth: 1

This page contains the documentation for the ``waveguide`` submodule, which functionally defines different kinds of waveguides, enables the creation of custom waveguides, and handles mesh generation for finite element analysis.

^^^^^^^^^^^^^^^^^^^
Waveguide class
^^^^^^^^^^^^^^^^^^^

.. autoclass:: wavesolve.waveguide.Waveguide
    :members:
    :exclude-members: plot_boundaries_recursive
    :member-order: groupwise
    
"""""""""""""""""""""""""""""""""""""""
classes that inherit from ``Waveguide``
"""""""""""""""""""""""""""""""""""""""

.. autoclass:: wavesolve.waveguide.CircularFiber
    :exclude-members: __new__

.. autoclass:: wavesolve.waveguide.EllipticalFiber
    :exclude-members: __new__

.. autoclass:: wavesolve.waveguide.PhotonicCrystalFiber
    :exclude-members: __new__

.. autoclass:: wavesolve.waveguide.PhotonicBandgapFiber
    :exclude-members: __new__

.. autoclass:: wavesolve.waveguide.FiberBundleLantern
    :exclude-members: __new__

.. autoclass:: wavesolve.waveguide.MCFPhotonicLantern
    :exclude-members: __new__

^^^^^^^^^^^^^^^^^^^
Prim2D class
^^^^^^^^^^^^^^^^^^^

.. autoclass:: wavesolve.waveguide.Prim2D
    :members:
    :exclude-members: get_nearest_bp_idx
    :member-order: groupwise

"""""""""""""""""""""""""""""""""""""""
classes that inherit from ``Prim2D``
"""""""""""""""""""""""""""""""""""""""

.. autoclass:: wavesolve.waveguide.Circle
    :exclude-members: __new__,__init__

    .. automethod:: make_points

.. autoclass:: wavesolve.waveguide.Ellipse
    :exclude-members: __new__,__init__

    .. automethod:: make_points

.. autoclass:: wavesolve.waveguide.Rectangle
    :exclude-members: __new__,__init__

    .. automethod:: make_points

.. autoclass:: wavesolve.waveguide.Prim2DUnion
    :exclude-members: __new__

.. autoclass:: wavesolve.waveguide.Prim2DArray
    :exclude-members: __new__