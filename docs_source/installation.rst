------------
Installation
------------

Make sure you have the following dependencies, e.g. using ``pip``: ::

    pip install numpy matplotlib gmsh pygmsh juliacall


``juliacall`` should automatically install a suitable version of Julia. [1]_


Then to install ``wavesolve``: ::

    pip install git+https://github.com/jw-lin/wavesolve.git@julia

Lastly, run the following (e.g. in IPython) to setup the Julia module ``FEsolver``: ::

    import wavesolve
    wavesolve.FEsolversetup()

This will set up the Julia dependencies ``SparseArrays``, ``PythonCall``, ``Arpack``, and ``Pardiso``, then precompile the Julia module.

To update the package: ::

    pip install --force-reinstall git+https://github.com/jw-lin/wavesolve.git@julia

and do the Julia module setup again.

If you are running into trouble with the Julia side, I recommend making sure ``juliacall`` is up to date, and/or emptying out the ``wavesolve`` (specifically the ``wavesolve/FEsolver/Manifest.toml`` file) before reinstalling.

.. rubric:: Footnotes

.. [1] For some reason, this tends to be an older version of Julia, at least on my machine.