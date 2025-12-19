------------
Installation
------------

To install ``wavesolve`` and its dependencies, use ``pip``: ::

    pip install git+https://github.com/jw-lin/wavesolve.git@julia

When importing ``wavesolve`` for the first time, a suitable Julia installation will be installed if necessary. Next run the following (e.g. in a Python shell) to setup the Julia module ``FEsolver`` and its dependencies: ::

    import wavesolve
    wavesolve.FEsolversetup()

To update the package: ::

    pip install --upgrade git+https://github.com/jw-lin/wavesolve.git@julia

and do the Julia module setup again. Now you should be able to run the snippets in :doc:`quickstart`.

Troubleshooting
---------------

If you are running into trouble with the Julia side, I recommend making sure ``juliacall`` is up to date, and/or emptying out the ``wavesolve`` folder (specifically the ``wavesolve/FEsolver/Manifest.toml`` file) before reinstalling.

For reference, the Python dependencies are ``numpy matplotlib gmsh pygmsh juliacall``, which are all ``pip``-installable.

The Julia dependencies are ``PythonCall LinearAlgebra SparseArrays Arpack``. To install these manually for a clean install of ``wavesolve`` you need to start a Julia prompt in the ``wavesolve`` directory, enter Julia's package manager (``]`` key), and then run the following: ::

    activate FEsolver
    add PythonCall LinearAlgebra SparseArrays Arpack
    precompile