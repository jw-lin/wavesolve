from juliacall import Pkg as jlPkg
import wavesolve,os

def FEsolversetup():
    path = os.path.dirname(wavesolve.__file__)
    jlPkg.activate(path+"/FEsolver")
    jlPkg.resolve()
    jlPkg.precompile()


