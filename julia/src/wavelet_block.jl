module WaveletBlock

using ..Utils
using ..WaveletOps
using ..Init
using Distributions
using Random
using StatsFuns: logistic

export update_cluster_wavelet_params_besov

function ensure_gamma_length!(wpar::WaveletParams, ncoeff::Int, M::Int)
    if length(wpar.gamma_ch) != M || any(length(g) != ncoeff for g in wpar.gamma_ch)
        wpar.gamma_ch = [Int.(rand(Bernoulli(0.2), ncoeff)) for _ in 1:M]
    end
end

function update_cluster_wavelet_params_besov(idx::Vector{Int}, precomp, M::Int, wpar::WaveletParams,
    sigma2_m::Vector{Float64}, tau_sigma::Float64;
    kappa_pi::Float64 = 0.6, c2::Float64 = 1.0, tau_pi::Float64 = 40.0,
    g_hyp = nothing,
    a_sig::Float64 = 2.5, b_sig::Float64 = 0.02,
    a_tau::Float64 = 2.0, b_tau::Float64 = 2.0,
    bias_coeff = nothing)

    if isempty(idx)
        return (wpar = wpar, beta_ch = wpar.gamma_ch, sigma2_m = sigma2_m, tau_sigma = tau_sigma, maps = nothing)
    end

    length(sigma2_m) == M || error("sigma2_m must be length M")
    stk = stack_D_from_precomp(precomp, idx, M; bias_coeff = bias_coeff)
    D = stk.D_arr
    maps = stk.maps
    ncoeff = size(D, 1)
    N = size(D, 2)
    lev_names = [String(k) for k in keys(maps[1].idx)]
    det_names = filter(name -> startswith(name, "d"), lev_names)
    s_name = filter(name -> startswith(name, "s"), lev_names)

    if isempty(wpar.pi_level)
        wpar.pi_level = Dict(name => 0.5 for name in det_names)
    end
    if isempty(wpar.g_level)
        wpar.g_level = Dict(name => 2.0 for name in det_names)
    end
    ensure_gamma_length!(wpar, ncoeff, M)

    # 1) gamma updates
    for m in 1:M
        Dm = view(D, :, :, m)
        gam = wpar.gamma_ch[m]
        for lev in det_names
            ids = maps[m].idx[Symbol(lev)]
            isempty(ids) && continue
            pi_j = wpar.pi_level[lev]
            g_j = wpar.g_level[lev]
            Dsub = Dm[ids, :]
            v_spike = sigma2_m[m]
            v_slab = (1 + g_j) * sigma2_m[m]
            ll_spike = -0.5 .* sum(log.(2π * v_spike) .+ (Dsub .^ 2) ./ v_spike; dims = 2)
            ll_slab = -0.5 .* sum(log.(2π * v_slab) .+ (Dsub .^ 2) ./ v_slab; dims = 2)
            logit_val = log(pi_j) .+ ll_slab .- (log(1 - pi_j) .+ ll_spike)
            p1 = logistic.(clamp.(logit_val, -35, 35))
            gam[ids] = Int.(rand.(Bernoulli.(vec(p1))))
        end
        if length(s_name) == 1
            ids_s = maps[m].idx[Symbol(s_name[1])]
            if !isempty(ids_s)
                gam[ids_s] .= 1
            end
        end
        wpar.gamma_ch[m] = gam
    end

    # 2) g_level updates
    for lev in det_names
        shape0 = isnothing(g_hyp) ? 2.0 : g_hyp[lev]["shape"]
        rate0 = isnothing(g_hyp) ? 2.0 : g_hyp[lev]["rate"]
        ss_over_sigma = 0.0
        n_sel_total = 0
        for m in 1:M
            ids = maps[m].idx[Symbol(lev)]
            isempty(ids) && continue
            Dm = view(D, :, :, m)
            sel = wpar.gamma_ch[m][ids] .== 1
            if any(sel)
                Did = Dm[ids[sel], :]
                ss_over_sigma += sum( Did .^ 2 ) / sigma2_m[m]
                n_sel_total += size(Did, 1)
            end
        end
        shape_post = shape0 + 0.5 * n_sel_total
        rate_post = rate0 + 0.5 * ss_over_sigma
        wpar.g_level[lev] = rand(InverseGamma(shape_post, 1 / rate_post))
    end

    # 3) pi_level updates
    for lev in det_names
        jnum = parse(Int, replace(lev, "d" => ""))
        m_j = clamp(kappa_pi * 2^(-c2 * jnum), 1e-6, 1 - 1e-6)
        a0 = tau_pi * m_j
        b0 = tau_pi * (1 - m_j)
        n1 = 0
        n0 = 0
        for m in 1:M
            ids = maps[m].idx[Symbol(lev)]
            isempty(ids) && continue
            gm = wpar.gamma_ch[m][ids]
            n1 += count(==(1), gm)
            n0 += count(==(0), gm)
        end
        wpar.pi_level[lev] = rand(Beta(a0 + n1, b0 + n0))
    end

    # 4) beta sampling
    beta_ch = [zeros(ncoeff) for _ in 1:M]
    for m in 1:M
        Dm = view(D, :, :, m)
        gam = wpar.gamma_ch[m]
        for lev in det_names
            ids = maps[m].idx[Symbol(lev)]
            isempty(ids) && continue
            g_j = wpar.g_level[lev]
            n = N
            Dbar = vec(mean(Dm[ids, :]; dims = 2))
            shrink = n / (n + 1 / g_j)
            mean_post = shrink .* Dbar
            var_post = sigma2_m[m] / (n + 1 / g_j)
            is_on = gam[ids] .== 1
            if any(is_on)
                dists = Normal.(mean_post[is_on], sqrt(var_post))
                beta_ch[m][ids[is_on]] = rand.(dists)
            end
        end
        if length(s_name) == 1
            ids_s = maps[m].idx[Symbol(s_name[1])]
            if !isempty(ids_s)
                Dbar_s = vec(mean(Dm[ids_s, :]; dims = 2))
                dists_s = Normal.(Dbar_s, sqrt(sigma2_m[m] / N))
                beta_ch[m][ids_s] = rand.(dists_s)
            end
        end
    end

    # 5) sigma2 updates & 6) tau_sigma
    sigma2_m_new = copy(sigma2_m)
    n_eff_m = ncoeff * N
    for m in 1:M
        Dm = view(D, :, :, m)
        resid = Dm .- beta_ch[m]
        ss_m = sum(resid .^ 2)
        shape_post = a_sig + 0.5 * n_eff_m
        rate_post = b_sig * tau_sigma + 0.5 * ss_m
        sigma2_m_new[m] = rand(InverseGamma(shape_post, 1 / rate_post))
    end
    a_post = a_tau + M * a_sig
    b_post = b_tau + b_sig * sum(1 ./ sigma2_m_new)
    tau_sigma_new = rand(Gamma(a_post, 1 / b_post))

    (wpar = wpar, beta_ch = beta_ch, sigma2_m = sigma2_m_new, tau_sigma = tau_sigma_new, maps = maps)
end

end # module
