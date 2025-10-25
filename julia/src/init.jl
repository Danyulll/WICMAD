module Init

using ..Utils
using ..Kernels
using ..WaveletOps
using ..ICMCache
using LinearAlgebra
using Random
using Distributions

export AcceptCount, AcceptanceTracker, WaveletParams, ClusterParams,
       draw_empty_acc, draw_new_cluster_params, ensure_complete_cache!

mutable struct AcceptCount
    a::Int
    n::Int
    AcceptCount() = new(0, 0)
end

mutable struct AcceptanceTracker
    kernel::Dict{Symbol,AcceptCount}
    L::AcceptCount
    eta::Vector{AcceptCount}
    tauB::AcceptCount
end

mutable struct WaveletParams
    lev_names::Vector{String}
    pi_level::Dict{String,Float64}
    g_level::Dict{String,Float64}
    gamma_ch::Vector{Vector{Int}}
end

mutable struct ClusterParams
    wpar::WaveletParams
    kern_idx::Int
    thetas::Vector{Dict{Symbol,Float64}}
    L::Matrix{Float64}
    eta::Vector{Float64}
    tau_B::Float64
    beta_ch::Vector{Vector{Float64}}
    sigma2::Vector{Float64}
    tau_sigma::Float64
    bias::Vector{Float64}
    cache::ICMCacheState
    mu_cached::Union{Nothing,Matrix{Float64}}
    mu_cached_iter::Int
    acc::AcceptanceTracker
    g_hyp::Union{Nothing,Dict{String,Dict{String,Float64}}}
end

function draw_empty_acc(M::Int, kernels::Vector{KernelConfig})
    pset = Set{Symbol}()
    for kc in kernels
        foreach(pn -> push!(pset, pn), kc.pnames)
    end
    kernel_counts = Dict{Symbol,AcceptCount}(pn => AcceptCount() for pn in pset)
    AcceptanceTracker(kernel_counts, AcceptCount(), [AcceptCount() for _ in 1:M], AcceptCount())
end

function draw_new_cluster_params(M::Int, P::Int, t, kernels::Vector{KernelConfig}; wf::String = "la8", J = nothing, boundary::String = "periodic")
    Jv = Utils.ensure_dyadic_J(P, J)
    zeros_mat = zeros(P, M)
    tmp = wt_forward_mat(zeros_mat; wf = wf, J = Jv, boundary = boundary)
    lev_names = [String(k) for k in keys(tmp[1].map.idx)]
    det_names = filter(name -> startswith(name, "d"), lev_names)
    ncoeff = length(tmp[1].coeff)
    beta_ch = [zeros(ncoeff) for _ in 1:M]
    pi_level = Dict(name => 0.5 for name in det_names)
    g_level = Dict(name => 2.0 for name in det_names)
    gamma_ch = [collect(rand(Binomial(1, 0.2), ncoeff)) for _ in 1:M]
    thetas = [kc.pstar() for kc in kernels]
    wpar = WaveletParams(lev_names, pi_level, g_level, gamma_ch)
    ClusterParams(
        wpar,
        rand(1:length(kernels)),
        thetas,
        Matrix{Float64}(I, M, M),
        fill(0.05, M),
        1.0,
        beta_ch,
        fill(0.05, M),
        1.0,
        zeros(M),
        ICMCacheState(),
        nothing,
        0,
        draw_empty_acc(M, kernels),
        nothing,
    )
end

function ensure_complete_cache!(cp::ClusterParams, kernels::Vector{KernelConfig}, t, M::Int)
    if length(cp.eta) != M || any(!isfinite.(cp.eta))
        cp.eta = fill(0.05, M)
    end
    if !isfinite(cp.tau_B) || cp.tau_B <= 0
        cp.tau_B = 1.0
    end
    if cp.kern_idx < 1 || cp.kern_idx > length(kernels)
        cp.kern_idx = 1
    end
    if length(cp.thetas) != length(kernels)
        cp.thetas = [kc.pstar() for kc in kernels]
    end
    if cp.thetas[cp.kern_idx] === nothing
        cp.thetas[cp.kern_idx] = kernels[cp.kern_idx].pstar()
    end
    if length(cp.bias) != M
        cp.bias = zeros(M)
    end
    cp.cache = build_icm_cache(t, kernels[cp.kern_idx], cp.thetas[cp.kern_idx], cp.L, cp.eta, cp.tau_B, cp.cache)
    cp
end

end # module
