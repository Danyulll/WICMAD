module Utils

using LinearAlgebra
using StatsBase
using StatsFuns: logit, logistic
using Distances
using Distributions

export scale_t01, dist_rows, nloc, ensure_dyadic_J, logit_safe, invlogit_safe,
       stick_to_pi, extend_sticks_until, update_v_given_z, pack_L, unpack_L,
       as_num_mat, normalize_t, eig_Kx, eig_Bshape, project_curve

function scale_t01(t::AbstractVector)
    mn = minimum(t); mx = maximum(t)
    rng = mx - mn
    rng <= 0 ? collect(t) : (collect(t) .- mn) ./ rng
end

function scale_t01(t::AbstractMatrix)
    T = Matrix{Float64}(t)
    for j in axes(T, 2)
        col = view(T, :, j)
        mn = minimum(col); mx = maximum(col)
        rng = mx - mn
        if rng > 0
            col .= (col .- mn) ./ rng
        end
    end
    T
end

scale_t01(t) = scale_t01(collect(t))

function dist_rows(t)
    if isa(t, AbstractVector)
        T = reshape(Float64.(t), :, 1)
    else
        T = Matrix{Float64}(t)
    end
    pairwise(Euclidean(), T; dims = 1)
end

nloc(t::AbstractVector) = length(t)
nloc(t::AbstractMatrix) = size(t, 1)

function ensure_dyadic_J(P::Integer, J::Union{Nothing, Integer, AbstractFloat})
    Jval = isnothing(J) ? log2(P) : Float64(J)
    Jint = round(Int, Jval)
    if abs(P - 2^Jint) > eps(Float64) * max(1, P)
        error("P must be 2^J (dyadic). Got P=$(P), J≈$(Jval) (rounded $(Jint) gives 2^J=$(2^Jint)).")
    end
    Jint
end

logit_safe(x) = logit(clamp(x, eps(Float64), 1 - eps(Float64)))
invlogit_safe(z) = logistic(z)

function stick_to_pi(v::AbstractVector)
    K = length(v)
    pi = zeros(Float64, K)
    tail = 1.0
    for k in 1:K
        pi[k] = v[k] * tail
        tail *= (1 - v[k])
    end
    pi
end

function extend_sticks_until(v::Vector{Float64}, alpha::Float64, threshold::Float64)
    tail = prod(1 .- v)
    while tail > threshold
        v_new = Base.rand(Beta(1, alpha))
        push!(v, v_new)
        tail *= (1 - v_new)
    end
    v
end

function update_v_given_z(v::Vector{Float64}, z::Vector{Int}, alpha::Float64)
    K = length(v)
    counts = counts(z, 1:K)
    tail_counts = reverse(cumsum(reverse(counts)))
    for k in 1:K
        a = 1 + counts[k]
        b = alpha + (k < K ? tail_counts[k + 1] : 0)
        v[k] = Base.rand(Beta(a, b))
    end
    v
end

function pack_L(L::AbstractMatrix)
    m, n = size(L)
    @assert m == n "L must be square"
    out = Vector{Float64}(undef, div(m * (m + 1), 2))
    idx = 1
    for j in 1:m
        for i in j:m
            out[idx] = L[i, j]
            idx += 1
        end
    end
    out
end

function unpack_L(theta::AbstractVector, m::Integer)
    expected = div(m * (m + 1), 2)
    length(theta) == expected || error("theta length $(length(theta)) incompatible with m=$m")
    L = zeros(Float64, m, m)
    idx = 1
    for j in 1:m
        for i in j:m
            L[i, j] = theta[idx]
            idx += 1
        end
    end
    for i in 1:m
        L[i, i] = abs(L[i, i]) + 1e-8
    end
    L
end

as_num_mat(A) = Matrix{Float64}(A)

function eig_Kx(Kx::AbstractMatrix)
    ee = eigen(Symmetric(Matrix{Float64}(Kx)))
    (V = ee.vectors, s = max.(ee.values, 0.0))
end

function eig_Bshape(L::AbstractMatrix, M::Int)
    Bshape = Matrix(L) * transpose(Matrix(L))
    trB = tr(Bshape)
    if trB > 0
        Bshape .= Bshape .* (M / trB)
    end
    ee = eigen(Symmetric(Bshape))
    (U = ee.vectors, lam = max.(ee.values, 0.0))
end

function project_curve(Yi::AbstractMatrix, mu::AbstractMatrix, V::AbstractMatrix)
    V === nothing && error("Eigenvectors V are nothing")
    transpose(V) * (Matrix{Float64}(Yi) - Matrix{Float64}(mu))
end

function normalize_t(t, P::Integer)
    if t === nothing
        error("t is nothing")
    end
    if isa(t, AbstractVector)
        length(t) == P || error("t vector has length $(length(t)) but P=$P")
        return collect(Float64.(t))
    end
    T = Matrix{Float64}(t)
    if size(T, 1) == P
        return T
    elseif size(T, 2) == P
        return transpose(T)
    else
        error("t has incompatible shape: $(size(T,1)) x $(size(T,2)); need P x d or length-P (P=$P)")
    end
end

end # module
