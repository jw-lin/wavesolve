from juliacall import Pkg as jlPkg
import wavesolve,os

def FEsolversetup():
    path = os.path.dirname(wavesolve.__file__)
    jlPkg.activate(path+"/FEsolver")
    jlPkg.add("PythonCall")
    jlPkg.add("SparseArrays")
    jlPkg.add("Arpack")
    jlPkg.precompile()
