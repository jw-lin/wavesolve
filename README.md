# wavesolve
`wavesolve` is a lightweight Python code to solve for the eigenmodes of waveguides. It can solve for both scalar and vector modes.
It uses the finite element method, and generates waveguide meshes through the `pygmsh` package. More details on the math behind `wavesolve` are included <a href="finite_element_method_notes.pdf">here</a>.

This particular branch moves some calculations into Julia; the Julia install is handled automatically in the below instructions. Meshing is still done in Python. Currently, this version solves the circular fiber waveguide example (scalar mode solving) from <a href="getting-started.ipynb">`getting-started.ipynb`</a> around 40% faster than base wavesolve. The vectorial solver is around $4\times$ faster judging from the doc examples.

Note that the first time the Julia-side code is called, things will run a few seconds slow. This is normal, and subsequent calls should be much faster.

## Installation

Make sure you have the following dependencies, e.g. using `pip`:

```
pip install numpy matplotlib gmsh pygmsh juliacall jupyter
```

`jupyter` is only needed to run the doc. The installation of `juliacall` will automatically trigger the installation of a suitable version of Julia (I think).[^1]

Then to install ``wavesolve``:

```
pip install git+https://github.com/jw-lin/wavesolve.git@julia
```

Lastly, run the following (e.g. in IPython) to setup the Julia module ``FEsolver``, which is intended to replace the old ``wavesolve.fe_solver``. Functions such as ``solve_waveguide`` are now under the ``FEsolver`` module.

```
import wavesolve
wavesolve.FEsolversetup()
```

This will set up the Julia dependencies ``SparseArrays``,``PythonCall``, ``Arpack``, and ``Pardiso``, then precompile the Julia module.

To update the package:
```
pip install --force-reinstall git+https://github.com/jw-lin/wavesolve.git@julia
```

and do the Julia module setup again.

If you are running into trouble with the Julia side, I recommend making sure ``juliacall`` is up to date, and/or emptying out the ``wavesolve`` (specifically the ``wavesolve/FEsolver/Manifest.toml`` file) before reinstalling.

## Documentation
See <a href="getting-started.ipynb">`getting-started.ipynb`</a> for an overview and some working examples. I may add a Sphinx web doc in the future.

## Acknowledgments
NSF grants 2109231, 2109232, 2308360, 2308361

[^1]: For some reason, ``juliacall`` doesn't want to install the latest version of Julia on my machine. But it seems fine.

