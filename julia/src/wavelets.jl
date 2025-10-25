module WaveletOps

using ..Utils
using Wavelets
using LinearAlgebra

export WaveletMap, WaveletCoefficients, wt_forward_1d, wt_inverse_1d, wt_forward_mat,
       precompute_wavelets, stack_D_from_precomp, compute_mu_from_beta

struct WaveletMap
    J::Int
    wf::String
    boundary::String
    P::Int
    idx::Dict{Symbol,UnitRange{Int}}
end

struct WaveletCoefficients
    coeff::Vector{Float64}
    map::WaveletMap
end

function wavelet_from_string(wf::String)
    sym = Symbol(lowercase(wf))
    if hasproperty(WT, sym)
        return wavelet(getproperty(WT, sym))
    else
        error("Unsupported wavelet family $wf")
    end
end

function extension_from_string(boundary::String)
    boundary == "periodic" && return Periodic()
    boundary == "reflection" && return Symmetric()
    error("Unsupported boundary $boundary; use \"periodic\" or \"reflection\".")
end

function wt_forward_1d(y::AbstractVector; wf::String = "la8", J::Union{Nothing,Int} = nothing, boundary::String = "periodic")
    P = length(y)
    Jv = Utils.ensure_dyadic_J(P, J)
    wave = wavelet_from_string(wf)
    ext = extension_from_string(boundary)
    wt = dwt(Float64.(y), wave, Jv; extension = ext)
    coeffs_vec = Float64[]
    idx = Dict{Symbol,UnitRange{Int}}()
    offset = 0
    for lev in 1:Jv
        d = detail(wt, lev)
        append!(coeffs_vec, d)
        idx[Symbol("d" * string(lev))] = offset + 1:offset + length(d)
        offset += length(d)
    end
    s = approx(wt)
    append!(coeffs_vec, s)
    idx[Symbol("s" * string(Jv))] = offset + 1:offset + length(s)
    WaveletCoefficients(coeffs_vec, WaveletMap(Jv, wf, boundary, P, idx))
end

function wt_inverse_1d(coeff_vec::AbstractVector, map::WaveletMap)
    wave = wavelet_from_string(map.wf)
    ext = extension_from_string(map.boundary)
    details = Vector{Vector{Float64}}(undef, map.J)
    for lev in 1:map.J
        key = Symbol("d" * string(lev))
        ids = map.idx[key]
        details[lev] = collect(coeff_vec[ids])
    end
    approx_coeffs = collect(coeff_vec[map.idx[Symbol("s" * string(map.J))]])
    wt = WaveletTransform(wave, approx_coeffs, details, ext)
    idwt(wt)
end

function wt_forward_mat(y_mat::AbstractMatrix; wf::String = "la8", J::Union{Nothing,Int} = nothing, boundary::String = "periodic")
    M = size(y_mat, 2)
    [wt_forward_1d(view(y_mat, :, m); wf = wf, J = J, boundary = boundary) for m in 1:M]
end

function precompute_wavelets(Y_list, wf::String, J, boundary::String)
    [wt_forward_mat(mat; wf = wf, J = J, boundary = boundary) for mat in Y_list]
end

function stack_D_from_precomp(precomp, idx::Vector{Int}, M::Int; bias_coeff = nothing)
    ncoeff = length(precomp[idx[1]][1].coeff)
    N = length(idx)
    D_arr = Array{Float64}(undef, ncoeff, N, M)
    for (jj, i) in enumerate(idx)
        for m in 1:M
            coeffs = precomp[i][m].coeff
            D_arr[:, jj, m] = coeffs
            if bias_coeff !== nothing
                bc = bias_coeff[m]
                if bc !== nothing && length(bc) == ncoeff
                    D_arr[:, jj, m] .-= bc
                end
            end
        end
    end
    maps = [precomp[idx[1]][m].map for m in 1:M]
    (; D_arr, maps)
end

function compute_mu_from_beta(beta_ch::Vector{Vector{Float64}}, wf::String, J::Int, boundary::String, P::Int)
    M = length(beta_ch)
    zeros_mat = zeros(P, M)
    tmpl = wt_forward_mat(zeros_mat; wf = wf, J = J, boundary = boundary)
    mu = zeros(P, M)
    for m in 1:M
        mu[:, m] = wt_inverse_1d(beta_ch[m], tmpl[m].map)
    end
    mu
end

end # module
