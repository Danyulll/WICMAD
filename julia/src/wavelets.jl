module WaveletOps

using ..Utils
using Wavelets
using Wavelets: WT
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

const _WAVELET_ALIASES = Dict(
    "haar" => :haar,
    "db1" => :haar,
)

function wavelet_from_string(wf::String, boundary::String = "periodic")
    sym = Symbol(lowercase(wf))
    sym = get(_WAVELET_ALIASES, String(sym), sym)
    if isdefined(WT, sym)
        ext = extension_from_string(boundary)
        return wavelet(getfield(WT, sym), WT.Filter, ext)
    else
        error("Unsupported wavelet family $wf")
    end
end

function extension_from_string(boundary::String)
    boundary == "periodic" && return WT.Periodic
    boundary == "reflection" && return WT.Symmetric
    error("Unsupported boundary $boundary; use \"periodic\" or \"reflection\".")
end

function wt_forward_1d(y::AbstractVector; wf::String = "la8", J::Union{Nothing,Int} = nothing, boundary::String = "periodic")
    P = length(y)
    Jv = Utils.ensure_dyadic_J(P, J)
    wave = wavelet_from_string(wf, boundary)
    wt = dwt(Float64.(y), wave, Jv)
    
    # The dwt function returns a flat vector of coefficients
    # We need to organize them into detail and approximation coefficients
    coeffs_vec = Float64[]
    idx = Dict{Symbol,UnitRange{Int}}()
    offset = 0
    
    # Extract detail coefficients for each level
    for lev in 1:Jv
        if lev == 1
            start_idx = 1
            end_idx = detailindex(P, 1, Jv) - 1
        else
            start_idx = detailindex(P, lev-1, Jv)
            end_idx = detailindex(P, lev, Jv) - 1
        end
        
        if start_idx <= end_idx
            d = wt[start_idx:end_idx]
            append!(coeffs_vec, d)
            idx[Symbol("d" * string(lev))] = offset + 1:offset + length(d)
            offset += length(d)
        end
    end
    
    # Extract approximation coefficients
    approx_start = detailindex(P, Jv, Jv)
    s = wt[approx_start:end]
    append!(coeffs_vec, s)
    idx[Symbol("s" * string(Jv))] = offset + 1:offset + length(s)
    
    WaveletCoefficients(coeffs_vec, WaveletMap(Jv, wf, boundary, P, idx))
end

function wt_inverse_1d(coeff_vec::AbstractVector, map::WaveletMap)
    wave = wavelet_from_string(map.wf, map.boundary)
    
    # Reconstruct the flat coefficient vector in the format expected by idwt
    wt_reconstructed = zeros(Float64, map.P)
    
    # Reconstruct detail coefficients
    for lev in 1:map.J
        key = Symbol("d" * string(lev))
        if haskey(map.idx, key)
            ids = map.idx[key]
            detail_coeffs = collect(coeff_vec[ids])
            
            # Place detail coefficients in the correct positions
            if lev == 1
                start_idx = 1
                end_idx = detailindex(map.P, 1, map.J) - 1
            else
                start_idx = detailindex(map.P, lev-1, map.J)
                end_idx = detailindex(map.P, lev, map.J) - 1
            end
            
            if start_idx <= end_idx && length(detail_coeffs) > 0
                wt_reconstructed[start_idx:end_idx] = detail_coeffs
            end
        end
    end
    
    # Reconstruct approximation coefficients
    approx_key = Symbol("s" * string(map.J))
    if haskey(map.idx, approx_key)
        ids = map.idx[approx_key]
        approx_coeffs = collect(coeff_vec[ids])
        approx_start = detailindex(map.P, map.J, map.J)
        if length(approx_coeffs) > 0
            wt_reconstructed[approx_start:end] = approx_coeffs
        end
    end
    
    # Perform inverse transform
    idwt(wt_reconstructed, wave, map.J)
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
