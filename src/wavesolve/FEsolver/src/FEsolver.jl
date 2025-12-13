module FEsolver
# construct matrices and solve eigenproblems in julia
using PythonCall
using SparseArrays
using Arpack

function compute_NN_dNdN(t::PyArray{Float64,2})
    return compute_NN_dNdN(pyconvert(Array{Float64,2},t))
end

function compute_NN_dNdN(t::Array{Float64,2}) :: Tuple{Vector{Float64},Vector{Float64}}
    NN = zeros(Float64,6,6)
    dNdN = zeros(Float64,6,6)

    x21 = t[1,2] - t[1,1]
    y21 = t[2,2] - t[2,1]
    x31 = t[1,3] - t[1,1]
    y31 = t[2,3] - t[2,1]
    x32 = t[1,3] - t[1,2]
    y32 = t[2,3] - t[2,2]
    _J = x21*y31 - x31*y21

    NN[1,1] = NN[2,2] = NN[3,3] = _J/60
    NN[4,4] = NN[5,5] = NN[6,6] = 4*_J/45
    NN[2,1] = NN[1,2] = NN[2,3] = NN[3,2] = NN[3,1] = NN[1,3] = -_J/360
    NN[1,5] = NN[5,1] = NN[2,6] = NN[6,2] = NN[3,4] = NN[4,3] = -_J/90
    NN[4,5] = NN[5,4] = NN[4,6] = NN[6,4] = NN[5,6] = NN[6,5] = 2*_J/45

    dNdN[1,1] = (y32*y32 + x32*x32)/(2*_J)
    dNdN[2,2] = (y31*y31 + x31*x31)/(2*_J)
    dNdN[3,3] = (y21*y21 + x21*x21)/(2*_J)
    dNdN[4,4] = 4/(3*_J) * (y32^2+y31*y21+x32^2+x31*x21)
    dNdN[5,5] = 4/(3*_J) * (y31^2-y21*y32+x31^2-x21*x32)
    dNdN[6,6] = 4/(3*_J) * (y21^2+y32*y31+x21^2+x32*x31)
    dNdN[1,2] = dNdN[2,1] = (y32*y31+x32*x31)/(6*_J)
    dNdN[2,3] = dNdN[3,2] = (y31*y21+x31*x21)/(6*_J)
    dNdN[3,1] = dNdN[1,3] = (-y21*y32-x21*x32)/(6*_J)
    dNdN[1,4] = dNdN[4,1] = dNdN[2,4] = dNdN[4,2]  = -2*(y32*y31+x32*x31)/(3*_J)
    dNdN[1,6] = dNdN[6,1] = dNdN[3,6] = dNdN[6,3] = 2*(y21*y32+x21*x32)/(3*_J)
    dNdN[3,5] = dNdN[5,3] = dNdN[2,5] = dNdN[5,2] = -2*(y31*y21+x31*x21)/(3*_J)
    dNdN[4,5] = dNdN[5,4] =  4/(3*_J)*(y21*y32+x21*x32)
    dNdN[5,6] = dNdN[6,5] = -4/(3*_J)*(y31*y32+x31*x32)
    dNdN[6,4] = dNdN[4,6] = -4/(3*_J)*(y31*y21+x31*x21)
    return vec(NN),vec(dNdN)
end

function construct_AB_order2_sparse(points:: Array{Float64,2},tris::Array{Int64,2}, IORs::Array{Float64,1}, k2::Float64) :: Tuple{SparseMatrixCSC{Float64, Int64}, SparseMatrixCSC{Float64, Int64}}
    N = size(points,2)
    Is,Js,Avals,Bvals = compute_AB_vals(points,tris,IORs,k2)
    A = sparse(Is,Js,Avals,N,N)
    B = sparse(Is,Js,Bvals,N,N)
    return A,B
end

function compute_AB_vals(points:: Array{Float64,2},tris::Array{Int64,2}, IORs::Array{Float64,1}, k2::Float64) :: Tuple{Vector{Int64},Vector{Int64},Vector{Float64},Vector{Float64}}
    Is = vec(repeat(tris,inner=[6,1]))
    Js = vec(repeat(tris,outer=[6,1]))

    Avals = Array{Float64}(undef,size(Is,1))
    Bvals = Array{Float64}(undef,size(Is,1))

    for i in axes(tris,2)
        n2 = IORs[i]
        tripoints = points[:,tris[:,i]]
        NN,dNdN = compute_NN_dNdN(tripoints)
        Avals[(i-1)*36+1:i*36] .= (NN*k2*n2 .- dNdN)
        Bvals[(i-1)*36+1:i*36] .= NN
    end
    return Is,Js,Avals,Bvals
end

function compute_AB_vals(points::PyArray{Float64,2},tris::PyArray{UInt64,2}, IORs::PyArray{Float64,1}, k2::Float64)
    return compute_AB_vals(pyconvert(Matrix{Float64},points),pyconvert(Matrix{Int64},tris),pyconvert(Vector{Float64},IORs),k2)
end

function construct_AB_order2_sparse(points:: PyArray{Float64,2},tris::PyArray{UInt64,2}, IORs::PyArray{Float64,1}, k2::Float64)
    return construct_AB_order2_sparse(pyconvert(Array{Float64,2},points),pyconvert(Array{Int64,2},tris),pyconvert(Array{Float64,1},IORs),k2)
end

function solve(A::SparseMatrixCSC,B::SparseMatrixCSC,est_eigval::Float64,Nmax::Int64)
    w,v = eigs(A, B, which=:LM,nev=Nmax,sigma=est_eigval,explicittransform=:none)
    return w,v
end

function solve_waveguide(points:: PyArray{Float64,2},tris::PyArray{UInt64,2}, IORs::PyArray{Float64,1}, k2::Float64,est_eigval::Float64,Nmax::Int64)
    A,B = construct_AB_order2_sparse(points,tris,IORs,k2)
    w,v = solve(A,B,est_eigval,Nmax)
    return w,v
end

end