# wavesolve

> **Update Dec 2025:** Check out the ``julia`` branch ([link](https://github.com/jw-lin/wavesolve/tree/julia)). Despite the name, it maintains the same Pythonic interface. Installation is slightly more involved and some functions are changed, but it's faster and has a basic web doc.

`wavesolve` is a lightweight Python code to solve for the eigenmodes of waveguides. It can solve for both scalar and vector modes.
It uses the finite element method, and generates waveguide meshes through the `pygmsh` package. More details on the math behind `wavesolve` are included <a href="finite_element_method_notes.pdf">here</a>.

## installation
Use pip: 

```
pip install git+https://github.com/jw-lin/wavesolve.git
```
To update the package:

```
pip install --force-reinstall git+https://github.com/jw-lin/wavesolve.git
```

Python dependencies: `numpy`,`scipy`,`matplotlib`,`numexpr`,`pygmsh`,`jupyter`,`pypardiso` (optional) \
Other dependencies: <a href="https://gmsh.info/">`Gmsh`</a> (required for `pygmsh`).

## documentation
See <a href="getting-started.ipynb">`getting-started.ipynb`</a> for an overview and some working examples.

## acknowledgments
NSF grants 2109231, 2109232, 2308360, 2308361



