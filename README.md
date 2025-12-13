# wavesolve
`wavesolve` is a lightweight Python code to solve for the eigenmodes of waveguides. It can solve for both scalar and vector modes.
It uses the finite element method, and generates waveguide meshes through the `pygmsh` package. More details on the math behind `wavesolve` are included <a href="finite_element_method_notes.pdf">here</a>.

This particular branch (work in progress) moves some calculations into Julia, and as such requires a Julia install.

## installation

Python dependencies: `numpy`,`scipy`,`matplotlib`,`numexpr`,`pygmsh`,`jupyter`, `juliacall`, `pypardiso` (optional, only useful for vector solving)

Other dependencies: <a href="https://gmsh.info/">`Gmsh`</a> (required for `pygmsh`).

Use pip: 

```
pip install git+https://github.com/jw-lin/wavesolve.git@julia
```

Then run the following (e.g. in IPython) to setup the Julia module ``FEsolver``, which is intended to replace the old ``wavesolve.fe_solver``. 

```
import wavesolve
wavesolve.FEsolversetup()
```

This will dowload the Julia dependencies ``SparseArrays``,``PythonCall``, and ``Arpack``, then precompile the Julia module.


To update the package:
```
pip install --force-reinstall git+https://github.com/jw-lin/wavesolve.git@julia
```

and do the Julia module setup again.

## documentation
See <a href="getting-started.ipynb">`getting-started.ipynb`</a> for an overview and some working examples. I may add a Sphinx web doc in the future.

## acknowledgments
NSF grants 2109231, 2109232, 2308360, 2308361



