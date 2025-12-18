# wavesolve
`wavesolve` is a lightweight Python code to solve for the eigenmodes of waveguides. It can solve for both scalar and vector modes.
It uses the finite element method, and generates waveguide meshes through the `pygmsh` package. More details on the math behind `wavesolve` are included <a href="finite_element_method_notes.pdf">here</a>.

This particular branch moves some calculations into Julia; the Julia install is handled automatically in the below instructions. Meshing is still done in Python. Currently, this version solves the circular fiber waveguide example from the doc around 40% faster than base wavesolve. The vectorial solver is around $4\times$ faster judging from the doc examples.

Note that the first time the Julia-side code is called, things will run a few seconds slow. This is normal, and subsequent calls should be much faster.

[Documentation](https://jw-lin.github.io/wavesolve/)

[How to install](https://jw-lin.github.io/wavesolve/installation.html)

Also, check out <a href="getting-started.ipynb">`getting-started.ipynb`</a> for the old notebook version of the doc.

## Acknowledgments
NSF grants 2109231, 2109232, 2308360, 2308361

[^1]: For some reason, ``juliacall`` doesn't want to install the latest version of Julia on my machine. But it seems fine.

