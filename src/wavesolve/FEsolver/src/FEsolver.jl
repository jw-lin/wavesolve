module FEsolver
# construct matrices and solve eigenproblems in julia
# also has some copy-pasted code from my cbeam package evaluating fields on FEMs
using PythonCall
using SparseArrays
using Arpack
using Pardiso

function solve_scalar(A::SparseMatrixCSC,B::SparseMatrixCSC,est_eigval::Float64,Nmax::Int64)
    w,v = eigs(A, B, which=:LM,nev=Nmax,sigma=est_eigval,explicittransform=:none)
    return w,v
end

function solve_vector(A::SparseMatrixCSC,B::SparseMatrixCSC,est_eigval::Float64,Nmax::Int64,solve_mode::String)
    if solve_mode == "transform"
        ps = MKLPardisoSolver()
        C = solve(ps,A .- est_eigval*B,Array(B))
        w,v = eigs(C,which=:SR,nev=Nmax)
        w = est_eigval .+ (1 ./ w)
    elseif solve_mode == "straight"
        w,v = eigs(A, B, which=:LM,nev=Nmax,sigma=est_eigval,explicittransform=:none)
    else
        w = nothing
        v = nothing
    end
    return w,v
end
#region scalar

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

function solve_waveguide(points:: PyArray{Float64,2},tris::PyArray{UInt64,2}, IORs::PyArray{Float64,1}, k2::Float64,est_eigval::Float64,Nmax::Int64)
    A,B = construct_AB_order2_sparse(points,tris,IORs,k2)
    w,v = solve_scalar(A,B,est_eigval,Nmax)
    return w,v
end

#endregion

#region vector

function solve_waveguide_vec(points:: PyArray{Float64,2},tris::PyArray{UInt64,2},edges::PyArray{UInt32,2},IORs::PyArray{Float64,1}, k2::Float64,Nedges::Int64, est_eigval::Float64,Nmax::Int64,solve_mode::String = "transform")
    A,B = construct_AB_vec(points,tris,edges,IORs,k2,Nedges)
    w,v = solve_vector(A,B,est_eigval,Nmax,solve_mode)
    return w,v
end

function construct_AB_vec(points:: PyArray{Float64,2},tris::PyArray{UInt64,2},edges::PyArray{UInt32,2},IORs::PyArray{Float64,1}, k2::Float64,Nedges::Int64)
    return construct_AB_vec( pyconvert(Array{Float64,2},points), pyconvert(Array{Int64,2},tris), pyconvert(Array{Int64,2},edges),pyconvert(Vector{Float64},IORs),k2::Float64,Nedges::Int64)
end

function construct_AB_vec(points:: Array{Float64,2},tris::Array{Int64,2},edges::Array{Int64,2},IORs::Array{Float64,1}, k2::Float64,Nedges::Int64)
    Ntt = Nedges
    Nzz = size(points,2)
    N = Ntt + Nzz
    Is_t = vec(repeat(edges,inner=[3,1]))
    Js_t = vec(repeat(edges,outer=[3,1]))

    Is_z = vec(repeat(tris,inner=[3,1]))
    Js_z = vec(repeat(tris,outer=[3,1]))

    Att = Array{Float64}(undef,size(Is_t,1))
    Btt = Array{Float64}(undef,size(Is_t,1))
    Btz = Array{Float64}(undef,size(Is_t,1))
    Bzt = Array{Float64}(undef,size(Is_t,1))
    Bzz = Array{Float64}(undef,size(Is_t,1))

    for i in axes(tris,2)
        n2 = IORs[i]
        tri_idx = tris[:,i]
        tripoints = points[:,tri_idx]
        pc = precompute(tripoints,tri_idx)
        NeNe = computeL_NeNe(pc)
        NN = computeL_NN(pc)
        dNdN = computeL_dNdN(pc)
        NedN = computeL_Ne_dN(pc)
        cdNcdN = computeL_curlNe_curlNe(pc)

        Att[(i-1)*9+1:i*9] .= n2*k2 * vec(NeNe) .- vec(cdNcdN)
        Btt[(i-1)*9+1:i*9] .= vec(NeNe)
        Btz[(i-1)*9+1:i*9] .= vec(NedN')
        Bzt[(i-1)*9+1:i*9] .= vec(NedN)
        Bzz[(i-1)*9+1:i*9] .= vec(dNdN) .- n2*k2*vec(NN)
    end

    # order is tt tz zt zz
    Is = [Is_t ; Is_t ; Is_z .+ Ntt ; Is_z .+ Ntt]
    Js = [Js_t ; Js_z .+ Ntt ; Js_t ; Js_z .+ Ntt]

    A = sparse(Is_t,Js_t,Att,N,N)
    B = sparse(Is,Js,[Btt ; Btz ; Bzt ; Bzz],N,N)
    return A,B
end
#endregion

#region finite element evaluation stuff

struct leaf
    bbox :: Array{Float64,1}
    idxs :: Vector{Int32}
end

struct idxtree
    left :: Union{idxtree,leaf,Nothing}
    right :: Union{idxtree,leaf,Nothing}
    bbox :: Union{Array{Float64,1},Nothing} # bbox: xmin,xmax,ymin,ymax
end

struct tritree # just need a way to store the array of triangle points
    _idxtree :: idxtree
    points :: Array{Float64,2} # array of all points in mesh
    connections :: Array{Int64,2}
    edges :: Array{Int64,2} # array of all edges, identified as a pair of indices for points
    tripoints :: Array{Float64,3}
end

function isleaf(a::Union{leaf,idxtree})
    return typeof(a) == leaf
end

function compute_SAH(cmins::AbstractVector{Float64},cmaxs::AbstractVector{Float64},csplit::Float64)
    total_interval = cmaxs[end]-cmins[1]
    area_A = (csplit-cmins[1])/total_interval
    area_B = 1 - area_A
    num_tris_A = sum(cmins .< csplit)
    num_tris_B = sum(cmaxs .> csplit)
    # assume node traversal costs 1 and point-tri intersect costs 4 
    #(empirical comparison between bbox test and triangle test)
    return 1 + area_A*num_tris_A*4 + area_B*num_tris_B*4
end

function minimize_SAH(cmins::AbstractVector{Float64},cmaxs::AbstractVector{Float64},num_buckets=8)
    start_minidx = Int(ceil(length(cmins)/2))
    N = size(cmins,1)
    if N <= 8
        return start_minidx # for small amounts of triangles SAH min isn't worth
    end
    default_cost = 4*N
    best_min_SAH = compute_SAH(cmins,cmaxs,cmins[start_minidx])

    splits = LinRange(cmins[1],cmins[end],num_buckets+2)[2:end-1]
    best_i = 4
    for i in 1:num_buckets
        s = splits[i]
        SAH = compute_SAH(cmins,cmaxs,s)
        if SAH <= best_min_SAH
            best_min_SAH = SAH
            best_i = i
        end
    end
    if best_min_SAH < default_cost
        return clamp(searchsortedfirst(cmins,splits[best_i])-1,1,N)
    else
        return nothing
    end
end

function computeBB(xmins,xmaxs,ymins,ymaxs)
    return [minimum(xmins),maximum(xmaxs),minimum(ymins),maximum(ymaxs)]
end

function computeBB(bounds::Matrix{Float64})::Vector{Float64}
    return [minimum(bounds[:,1]),maximum(bounds[:,2]),minimum(bounds[:,3]),maximum(bounds[:,4])]
end

function construct_recursive(idxs::Array{Int64,1},bounds::Array{Float64,2},level::Int64,min_leaf_size::Int64)
    BB = computeBB(bounds)
    if size(idxs,1) <= min_leaf_size
        return leaf(BB,idxs)
    else
        ax = level%2+1
        sort_idx = sortperm(bounds[:,(ax-1)*2+1])
        cmins = bounds[sort_idx,2*ax-1]
        cmaxs = bounds[sort_idx,2*ax]

        split_idx = minimize_SAH(cmins,cmaxs)
        if isnothing(split_idx)
            return leaf(BB,idxs)
        end

        left_idxs = sort_idx[1:split_idx]
        right_idxs = sort_idx[split_idx+1:end]

        if !isempty(left_idxs)                
            left_bounds = bounds[left_idxs,:]
            left_tree = construct_recursive(idxs[left_idxs],left_bounds,level+1,min_leaf_size)
        else
            left_tree = nothing
        end

        if !isempty(right_idxs)
            right_bounds = bounds[right_idxs,:]
            right_tree = construct_recursive(idxs[right_idxs],right_bounds,level+1,min_leaf_size)
        else
            right_tree = nothing
        end
        return idxtree(left_tree,right_tree,BB)
    end
end

function construct_tritree(points::PyArray{Float64,2},connections::PyArray{T,2} where T<:Integer,edges::PyArray{T,2} where T<:Integer,min_leaf_size=4)
    points = pyconvert(Array{Float64,2},points)
    connections = pyconvert(Array{UInt64,2},connections)

    tripoints = Array{Float64,3}(undef,size(connections,1),size(connections,2),2)
    for i in axes(connections,1)
        for j in axes(connections,2)
            tripoints[i,j,1] = points[connections[i,j],1]
            tripoints[i,j,2] = points[connections[i,j],2]
        end
    end

    xmins = minimum(tripoints[:,:,1],dims=2)
    xmaxs = maximum(tripoints[:,:,1],dims=2)
    ymins = minimum(tripoints[:,:,2],dims=2)
    ymaxs = maximum(tripoints[:,:,2],dims=2)
    bounds = hcat(xmins,xmaxs,ymins,ymaxs)

    idxs = Array(1:size(tripoints,1))
    return tritree(construct_recursive(idxs,bounds,0,min_leaf_size),points,connections,edges,tripoints)
end 

function inside(point::AbstractVector{Float64},bbox::AbstractVector{Float64})
    return (bbox[1]<=point[1]<=bbox[2]) & (bbox[3]<=point[2]<=bbox[4])
end

function det(u::AbstractVector{Float64},v::AbstractVector{Float64})
    return u[1]*v[2] - u[2]*v[1]
end

function inside(point::AbstractVector{Float64},tri::Array{Float64,2},_eps=1e-12)::Bool
    x,y = point
    dot1 = (tri[2,2]-tri[1,2])*(x-tri[1,1]) + (tri[1,1]-tri[2,1])*(y-tri[1,2])
    if dot1 > _eps
        return false
    end
    dot2 = (tri[3,2]-tri[2,2])*(x-tri[2,1]) + (tri[2,1]-tri[3,1])*(y-tri[2,2])
    if dot2 > _eps
        return false
    end
    dot3 = (tri[1,2]-tri[3,2])*(x-tri[3,1]) + (tri[3,1]-tri[1,1])*(y-tri[3,2])
    if dot3 > _eps
        return false
    end
    return true
end

function query_recursive(point::Union{AbstractVector{Float64},PyArray{Float64,1}},tripoints::Array{Float64,3},_idxtree::Union{idxtree,leaf,Nothing}) :: Int64
    if typeof(point) <: PyArray
        point = pyconvert(Vector{Float64},point)
    end
    if isnothing(_idxtree)
        return 0
    elseif !inside(point,_idxtree.bbox)
        return 0
    elseif isleaf(_idxtree)
        for idx in _idxtree.idxs
            tri = tripoints[idx,:,:]
            if inside(point,tri)
                return idx
            end
        end
        return 0
    else
        left_query = query_recursive(point,tripoints,_idxtree.left)
        if left_query != 0
            return left_query
        end
        right_query = query_recursive(point,tripoints,_idxtree.right)
        if right_query != 0
            return right_query
        else
            return 0
        end
    end
end

function query(point::Union{AbstractVector{Float64},PyArray{Float64,1}},_tritree::tritree) :: Int64
    if typeof(point) <: PyArray
        point = pyconvert(Vector{Float64},point)
    end
    return query_recursive(point,_tritree.tripoints,_tritree._idxtree)
end

### interpolation stuff ##

function affine_transform_matrix(vertices::AbstractArray{Float64,2} )
    x21 = vertices[2,1] - vertices[1,1]
    y21 = vertices[2,2] - vertices[1,2]
    x31 = vertices[3,1] - vertices[1,1]
    y31 = vertices[3,2] - vertices[1,2]
    _J = x21*y31-x31*y21
    M = [y31 -x31 ; -y21 x21] ./ _J
    return M
end

function affine_transform_matrixT_in_place(vertices::Union{Matrix{Float64},SubArray{Float64,2}},M::Matrix{Float64} )
    x21 = vertices[2,1] - vertices[1,1]
    y21 = vertices[2,2] - vertices[1,2]
    x31 = vertices[3,1] - vertices[1,1]
    y31 = vertices[3,2] - vertices[1,2]
    _J = x21*y31-x31*y21
    M[1,1] = y31/_J
    M[2,1] = -x31/_J
    M[1,2] = -y21/_J
    M[2,2] = x21/_J
    return M
end

function affine_transform_matrix_inv(vertices::Array{Float64,2})
    x21 = vertices[2,1] - vertices[1,1]
    y21 = vertices[2,2] - vertices[1,2]
    x31 = vertices[3,1] - vertices[1,1]
    y31 = vertices[3,2] - vertices[1,2]
    return [x21 x31 ; y21 y31]
end

function apply_affine_transform(vertices::AbstractArray{Float64,2},xy::AbstractVector{Float64})
    M = affine_transform_matrix(vertices)
    return M * (xy .- vertices[1,:])
end

function jac(vertices::Array{Float64,2})
    x21 = vertices[2,1] - vertices[1,1]
    y21 = vertices[2,2] - vertices[1,2]
    x31 = vertices[3,1] - vertices[1,1]
    y31 = vertices[3,2] - vertices[1,2]
    return x21*y31-x31*y21
end

function get_interp_weights(new_points::PyArray{Float64,2},_tritree::tritree)
    _old_points = _tritree.tripoints
    _new_points = pyconvert(Array{Float64,2},new_points)
    N = size(_new_points,1)
    _weights = Array{Float64}(undef,N,6)
    _triidxs = Array{Int64}(undef,N)

    for i in 1:N
        new_point = _new_points[i,1:2]
        _triidx = query(new_point,_tritree)
        if _triidx != 0
            _triidxs[i] = _triidx
            triverts = _old_points[_triidx,:,:]
            u,v = apply_affine_transform(triverts,new_point)

            _weights[i,1] = N1(u,v)
            _weights[i,2] = N2(u)
            _weights[i,3] = N3(v)
            _weights[i,4] = N4(u,v)
            _weights[i,5] = N5(u,v)
            _weights[i,6] = N6(u,v)
        else
            _weights[i,:] .= 0.
            _triidxs[i] = 0
        end
    end
    return (_triidxs,_weights)
end

function evaluate_vec(point::Union{AbstractVector{Float64},PyArray{Float64,1}},field::Union{PyArray{T,1},Vector{T}},_tritree::tritree) :: Vector{T} where T<:Union{Float64,ComplexF64}
    """pointwise evaulation of a vectorial field"""
    dtype = eltype(field)
    if typeof(point) <: PyArray
        point = pyconvert(Vector{Float64},point)
    end
    if typeof(field) <: PyArray
        field = pyconvert(Vector{dtype},field)
    end
    _triidx = query(point,_tritree)
    _triidxs = _tritree.connections[_triidx,:]
    _edgeidxs = _tritree.edges[_triidx,:]
    val = [0. , 0.]
    if _triidx != 0
        triverts = _tritree.tripoints[_triidx,:,:]
        val .+= LNe0(point, triverts, _triidxs) * field[_edgeidxs[1]]
        val .+= LNe1(point, triverts, _triidxs) * field[_edgeidxs[2]]
        val .+= LNe2(point, triverts, _triidxs) * field[_edgeidxs[3]]
    end
    return val
end

function evaluate(point::Union{AbstractVector{Float64},PyArray{Float64,1}},field::Union{PyArray{T,1},Vector{T}},_tritree::tritree ; order::Int64 = 2) :: T where T<:Union{Float64,ComplexF64}
    """pointwise evaulation of a scalar field"""
    dtype = eltype(field)
    if typeof(point) <: PyArray
        point = pyconvert(Vector{Float64},point)
    end
    if typeof(field) <: PyArray
        field = pyconvert(Vector{dtype},field)
    end
    _triidx = query(point,_tritree)
    val = 0.
    if _triidx != 0
        @views triverts = _tritree.tripoints[_triidx,:,:]
        u,v = apply_affine_transform(triverts,point)
        if order==2
            val += N1(u,v) * field[_tritree.connections[_triidx,1]]
            val += N2(u) * field[_tritree.connections[_triidx,2]]
            val += N3(v) * field[_tritree.connections[_triidx,3]]
            val += N4(u,v) * field[_tritree.connections[_triidx,4]]
            val += N5(u,v) * field[_tritree.connections[_triidx,5]]
            val += N6(u,v) * field[_tritree.connections[_triidx,6]]
        elseif order==1
            val += LN1(u,v) * field[_tritree.connections[_triidx,1]]
            val += LN2(u,v) * field[_tritree.connections[_triidx,2]]
            val += LN3(u,v) * field[_tritree.connections[_triidx,3]]
        end
    end
    return val
end

function evaluate(point::Union{PyArray{Float64,2},Matrix{Float64}},field::Union{PyArray{T,1},Vector{T}},_tritree::tritree ; order::Int64 = 2) :: Vector{T} where T<:Union{Float64,ComplexF64}
    """linewise evaluation of a scalar field"""
    dtype = eltype(field)
    if typeof(point) <: PyArray
        point = pyconvert(Matrix{Float64},point)
    end
    if typeof(field) <: PyArray
        field = pyconvert(Vector{dtype},field)
    end
    out = Vector{dtype}(undef,size(point,1))

    for i in axes(point,1)
        @views _point = point[i,:]
        out[i] = evaluate(_point,field,_tritree,order=order)
    end
    return out
end

function evaluate_vec(point::Union{PyArray{Float64,2},Matrix{Float64}},field::Union{PyArray{T,1},Vector{T}},_tritree::tritree) :: Array{T,2} where T<:Union{Float64,ComplexF64}
    """linewise evaluation of a vector field"""
    dtype = eltype(field)
    if typeof(point) <: PyArray
        point = pyconvert(Matrix{Float64},point)
    end
    if typeof(field) <: PyArray
        field = pyconvert(Vector{dtype},field)
    end
    out = Array{dtype}(undef,size(point,1),2)

    for i in axes(point,1)
        @views _point = point[i,:] # technically further in memory i think
        out[i,:] .= evaluate_vec(_point,field,_tritree)
    end
    return out
end

function evaluate(pointsx::PyArray{Float64,1},pointsy::PyArray{Float64,1},field::PyArray{T,1} where T<:Union{Float64,ComplexF64},_tritree::tritree ; order::Int64 = 2)
    """gridwise evaluation of a scalar field"""
    dtype = eltype(field)
    pointsx = pyconvert(Vector{Float64},pointsx)
    pointsy = pyconvert(Vector{Float64},pointsy)
    field = pyconvert(Vector{dtype},field)
    out = Array{dtype,2}(undef,size(pointsx,1),size(pointsy,1))

    for i in eachindex(pointsx)
        for j in eachindex(pointsy)
            point = [pointsx[i],pointsy[j]]
            out[i,j] = evaluate(point,field,_tritree,order=order)
        end
    end
    return out
end

function evaluate_vec(pointsx::PyArray{Float64,1},pointsy::PyArray{Float64,1},field::PyArray{T,1} where T<:Union{Float64,ComplexF64},_tritree::tritree)
    """gridwise evaluation of a vector field"""
    dtype = eltype(field)
    pointsx = pyconvert(Vector{Float64},pointsx)
    pointsy = pyconvert(Vector{Float64},pointsy)
    field = pyconvert(Vector{dtype},field)
    out = Array{dtype,3}(undef,size(pointsx,1),size(pointsy,1),2)

    for i in eachindex(pointsx)
        for j in eachindex(pointsy)
            point = [pointsx[i],pointsy[j]]
            out[i,j,:] .= evaluate_vec(point,field,_tritree)
        end
    end
    return out
end

function evaluate_func(field::PyArray{T,1} where T<:Union{Float64,ComplexF64},_tritree::tritree)
    """ convert a FE field into a function of point [x,y] """
    dtype = eltype(field)
    field = pyconvert(Vector{dtype},field)

    function _inner_(point::Union{AbstractVector{Float64},PyArray{Float64,1}})
        if typeof(point) <: PyArray
            point = pyconvert(Vector{Float64},point)
        end
        return evaluate(point,field,_tritree)
    end
    return _inner_
end

function evaluate(f::Function,xa::PyArray{Float64,1},ya::PyArray{Float64,1})
    xa = pyconvert(Vector{Float64},xa)
    ya = pyconvert(Vector{Float64},ya)
    out = Array{Float64}(undef,size(xa,1),size(ya,1))
    for i in eachindex(xa)
        for j in eachindex(ya)
            x = xa[i]
            y = ya[j]
            out[i,j] = f([x,y])
        end
    end
    return out
end

#endregion

#region shapefuncs - scalar
function N1(u,v)
    2 * (1 - u - v) * (0.5 - u - v)
end
function N2(u)
    2 * u * (u - 0.5)
end
function N3(v)
    2 * v * (v - 0.5)
end
function N4(u,v)
    4 * u * (1 - u - v)
end
function N5(u,v)
    4 * u * v
end
function N6(u,v)
    4 * v * (1 - u - v)
end
function duN1(u,v)
    4 * (u + v) - 3
end
function dvN1(u,v)
    4 * (u + v) - 3
end
function duN2(u,v)
    4 * u - 1
end
function dvN2(u,v)
    0.
end
function duN3(u,v)
    0.
end
function dvN3(u,v)
    4 * v - 1
end
function duN4(u,v)
    -8*u - 4*v + 4
end
function dvN4(u,v)
    -4 * u
end
function duN5(u,v)
    4 * v
end
function dvN5(u,v)
    4 * u
end
function duN6(u,v)
    -4 * v
end
function dvN6(u,v)
    -4*u - 8*v + 4
end
#endregion

#region shapefuncs-scalar linear
function LN1(u,v)
    1 - u - v
end
function LN2(u,v)
    u
end
function LN3(u,v)
    v
end

#region shapefuncs-vector

function precompute(tri::Array{Float64,2},tri_idx::Vector{Int}) :: NTuple{10,Float64}
    # compute a bunch of stuff for shapefuncs
    x21 = tri[1,2] - tri[1,1]
    y21 = tri[2,2] - tri[2,1]
    x31 = tri[1,3] - tri[1,1]
    y31 = tri[2,3] - tri[2,1]
    x32 = tri[1,3] - tri[1,2]
    y32 = tri[2,3] - tri[2,2]

    s1 = (tri_idx[1]<tri_idx[2] ? 1 : -1)
    s2 = (tri_idx[2]<tri_idx[3] ? 1 : -1)
    s3 = (tri_idx[3]<tri_idx[1] ? 1 : -1)

    l12 = sqrt(x21*x21+y21*y21)*s1 
    l23 = sqrt(x32*x32+y32*y32)*s2
    l31 = sqrt(x31*x31+y31*y31)*s3
    _J = x21*y31 - x31*y21

    return (x21,y21,x31,y31,x31,y31,l12,l23,l31,_J)
end

function precompute_row_ordered(tri::Union{Array{Float64,2},PyArray{Float64,2}},tri_idx::Union{Vector{Int},PyArray{UInt64,1}}) :: NTuple{10,Float64}
    x21 = tri[2,1] - tri[1,1]
    y21 = tri[2,2] - tri[1,2]
    x31 = tri[3,1] - tri[1,1]
    y31 = tri[3,2] - tri[1,2]
    x32 = tri[3,1] - tri[2,1]
    y32 = tri[3,2] - tri[2,2]

    s1 = (tri_idx[1]<tri_idx[2] ? 1 : -1)
    s2 = (tri_idx[2]<tri_idx[3] ? 1 : -1)
    s3 = (tri_idx[3]<tri_idx[1] ? 1 : -1)

    l12 = sqrt(x21*x21+y21*y21)*s1 
    l23 = sqrt(x32*x32+y32*y32)*s2
    l31 = sqrt(x31*x31+y31*y31)*s3
    _J = x21*y31 - x31*y21

    return (x21,y21,x31,y31,x31,y31,l12,l23,l31,_J)
end

function LNe0(p,tri,tri_idx)
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precompute_row_ordered(tri,tri_idx)
    x1,y1 = tri[1,1],tri[1,2]
    x,y = p[1],p[2]
    xv = (-y + y31 + y1)/_J*l12
    yv = (x - x31 - x1)/_J*l12
    return [xv , yv]
end

function LNe1(p,tri,tri_idx)
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precompute_row_ordered(tri,tri_idx)
    x1,y1 = tri[1,1],tri[1,2]
    x,y = p[1],p[2]
    xv = (-y + y1)/_J*l23
    yv = (x - x1)/_J*l23
    return [xv , yv]
end

function LNe2(p,tri,tri_idx)
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precompute_row_ordered(tri,tri_idx)
    x1,y1 = tri[1,1],tri[1,2]
    x,y = p[1],p[2]
    xv = (-y + y21 + y1)/_J*l31
    yv = (x - x21 - x1)/_J*l31
    return [xv , yv]
end

function computeL_NeNe(precomp::NTuple{10,Float64})
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precomp
    out = zeros(Float64,3,3)
    out[1,1] = l12^2 * (l12^2-3*x21*x31+3*l31^2-3*y21*y31)/(12*_J)
    out[2,2] = l23^2 * (l12^2 + x21*x31 + l31^2 + y21*y31)/(12*_J)
    out[3,3] = l31^2 * (3*l12^2-3*x21*x31+l31^2-3*y21*y31)/(12*_J)
    out[1,2] = out[2,1] = l12*l23 * (l12^2-x21*x31-l31^2-y21*y31)/(12*_J)
    out[2,3] = out[3,2] = l23*l31 * (-l12^2-x21*x31+l31^2-y21*y31)/(12*_J)
    out[3,1] = out[1,3] = l31*l12 * (-l12^2+3*x21*x31-l31^2+3*y21*y31)/(12*_J)
    return out
end

function computeL_curlNe_curlNe(precomp::NTuple{10,Float64})
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precomp
    ls = [l12 , l23 , l31]
    return 2/_J * ls * ls'
end

function computeL_Ne_dN(precomp::NTuple{10,Float64})
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precomp
    out = zeros(Float64,3,3)
    out[1,1] = l12*(-l12^2 + 3*x21*x31-2*l31^2+3*y21*y31)/(6*_J)
    out[2,2] = l23*(-x21*x31-l31^2-y21*y31)/(6*_J)
    out[3,3] = l31*(-2*l12^2+x21*x31+y21*y31)/(6*_J)
    out[1,2] = l12*(-x21*x31+2*l31^2-y21*y31)/(6*_J)
    out[2,1] = l23*(-l12^2+l31^2)/(6*_J)
    out[2,3] = l23*(l12^2+x21*x31+y21*y31)/(6*_J)
    out[3,2] = l31*(2*x21*x31+2*y21*y31-l31^2)/(6*_J)
    out[3,1] = l31*(2*l12^2-3*x21*x31+l31^2-3*y21*y31)/(6*_J)
    out[1,3] = l12*(l12^2-2*x21*x31-2*y21*y31)/(6*_J)
    return out
end

function computeL_NN(precomp::NTuple{10,Float64})
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precomp
    out = fill(_J/24,3,3)
    out[1,1] = out[2,2] = out[3,3] = _J/12
    return out
end

function computeL_dNdN(precomp::NTuple{10,Float64})
    x21,y21,x31,y31,x31,y31,l12,l23,l31,_J = precomp
    out = zeros(Float64,3,3)
    out[1,1] = (l12^2+l31^2-2*x21*x31-2*y21*y31)/(2*_J)
    out[2,2] = l31^2/(2*_J)
    out[3,3] = l12^2/(2*_J)
    out[1,2] = out[2,1] = (-l31^2+x21*x31+y21*y31)/(2*_J)
    out[2,3] = out[3,2] = (-x21*x31-y21*y31)/(2*_J)
    out[3,1] = out[1,3] = (-l12^2+x21*x31+y21*y31)/(2*_J)
    return out
end

#endregion

end # module FEsolver