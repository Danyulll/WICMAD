module Kernels

using ..Utils
using LinearAlgebra
using StatsFuns: logit
using Distributions

export KernelConfig, make_kernels

struct KernelConfig
    name::String
    fun::Function
    pnames::Vector{Symbol}
    prior::Function
    pstar::Function
    prop_sd::Dict{Symbol,Float64}
end

function k_sqexp(t, l_scale)
    D2 = Utils.dist_rows(t).^2
    @. exp(-0.5 * D2 / (l_scale^2))
end

function k_mat32(t, l_scale)
    D = Utils.dist_rows(t)
    r = D ./ l_scale
    a = sqrt(3.0) .* r
    @. (1 + a) * exp(-a)
end

function k_mat52(t, l_scale)
    D = Utils.dist_rows(t)
    r = D ./ l_scale
    a = sqrt(5.0) .* r
    @. (1 + a + 5 * r^2 / 3) * exp(-a)
end

function k_periodic(t, l_scale, period)
    D = Utils.dist_rows(t)
    @. exp(-2 * sinpi(D / period)^2 / (l_scale^2))
end

k_bias(t, s0) = (s0^2) .* ones(Float64, Utils.nloc(t), Utils.nloc(t))

function make_bias_variant(kcfg::KernelConfig)
    fun = function(t, par)
        base = kcfg.fun(t, par)
        base + k_bias(t, par[:s0])
    end
    prior = function(par)
        kcfg.prior(par) + logpdf(Normal(-3, 0.75), log(max(par[:s0], 1e-9)))
    end
    pstar = function()
        base = kcfg.pstar()
        base[:s0] = exp(Base.rand(Normal(-3, 0.75)))
        base
    end
    prop_sd = copy(kcfg.prop_sd)
    prop_sd[:s0] = 0.25
    KernelConfig(kcfg.name * "+Bias", fun, vcat(kcfg.pnames, :s0), prior, pstar, prop_sd)
end

function make_kernels(; add_bias_variants::Bool = true)
    base = KernelConfig[
        KernelConfig(
            "SE",
            (t, par) -> k_sqexp(t, par[:l_scale]),
            [:l_scale],
            par -> logpdf(Gamma(2, 1 / 2), par[:l_scale]),
            () -> Dict(:l_scale => Base.rand(Gamma(2, 1 / 2))),
            Dict(:l_scale => 0.20),
        ),
        KernelConfig(
            "Mat32",
            (t, par) -> k_mat32(t, par[:l_scale]),
            [:l_scale],
            par -> logpdf(Gamma(2, 1 / 2), par[:l_scale]),
            () -> Dict(:l_scale => Base.rand(Gamma(2, 1 / 2))),
            Dict(:l_scale => 0.20),
        ),
        KernelConfig(
            "Mat52",
            (t, par) -> k_mat52(t, par[:l_scale]),
            [:l_scale],
            par -> logpdf(Gamma(2, 1 / 2), par[:l_scale]),
            () -> Dict(:l_scale => Base.rand(Gamma(2, 1 / 2))),
            Dict(:l_scale => 0.20),
        ),
        KernelConfig(
            "Periodic",
            (t, par) -> k_periodic(t, par[:l_scale], par[:period]),
            [:l_scale, :period],
            par -> logpdf(Gamma(3, 1 / 2), par[:l_scale]) + logpdf(Beta(5, 5), par[:period]),
            () -> Dict(:l_scale => Base.rand(Gamma(3, 1 / 2)), :period => Base.rand(Beta(5, 5))),
            Dict(:l_scale => 0.20, :period => 0.20),
        ),
    ]

    if !add_bias_variants
        return base
    end

    vcat(base, [make_bias_variant(k) for k in base])
end

end # module
