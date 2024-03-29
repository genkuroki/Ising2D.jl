"""
Ising2D is a simple Julia module of 2D Ising model.

Benchmark 1:

```julia
using Ising2D, BenchmarkTools
s = rand_ising2d()
print("VERSION = \$VERSION:")
@btime ising2d!(\$s; algorithm=IfElse());
```

Benchmark 2:

```julia
using Ising2D, BenchmarkTools
s = rand_ising2d()
print("VERSION = \$VERSION:")
@btime ising2d!(\$s; algorithm=MultiFor());
```
"""
module Ising2D

using Plots
using Plots.PlotMeasures
using Random
using Random: default_rng
using ProgressMeter

export
    default_rng, 
    β_ising2d, rand_ising2d, ones_ising2d,
    IfElse, MultiFor, ising2d!,
    plot_ising2d, gif_ising2d,
    mcmc_ising2d!, energy_density_ising2d, magnetization_ising2d,
    plot_mcmc_ising2d, histogram_mcmc_ising2d, gif_mcmc_ising2d

"""
    β_ising2d = log(1 + √2)/2

is the critical inverse temperature of 2D Ising model on the square lattice.
"""
const β_ising2d = log(1 + √2)/2

"""
    rand_ising2d(rng::AbstractRNG, m=100, n=m)
    rand_ising2d(m=100, n=m)

generate the random state of 2D Ising model.
"""
rand_ising2d(rng::AbstractRNG, m=100, n=m) = rand(rng, Int8[-1, 1], m, n)
rand_ising2d(m=100, n=m) = rand(Int8[-1, 1], m, n)

"""
    ones_ising2d(m=100, n=m)

generates the all-one state of 2D Ising model.
"""
ones_ising2d(m=100, n=m) = ones(Int8, m, n)

"""
`IfElse()` = Algorithm using ifelse
"""
struct IfElse end

"""
`MultiFor()` = Algorithm using multiple for-loops
"""
struct MultiFor end

"""
`default_algorithm() = MultiFor()`
"""
default_algorithm() = MultiFor()

"""
    ising2d!(s=rand_ising2d(), β=β_ising2d, niters=10^3, rng=default_rng();
        algorithm=default_algorithm())

updates the 2-dimensional array `s` with values ±1 randomly by 2D Ising model rule.

Example 1:

```julia
using Ising2D
@time s = ising2d!(rand_ising2d(), β_ising2d, 10^5)
plot_ising2d(s)
```

Example 2:

```julia
using Ising2D
@time s = ising2d!(rand_ising2d(), β_ising2d, 10^5; algorithm=IfElse())
plot_ising2d(s)
```
"""
function ising2d!(s=rand_ising2d(), β=β_ising2d, niters=10^3, rng=default_rng();
        algorithm=default_algorithm())
    ising2d!(algorithm, s, β, niters, rng)
end

"""
    ising2d!(::IfElse, s=rand_ising2d(), β=β_ising2d, niters=10^3, rng=default_rng())

uses the algorithm using ifelse in only one for-loop.

Example:

```julia
s = ones_ising2d()
ising2d!(IfElse(), s)
plot_ising2d(s)
```
"""
function ising2d!(::IfElse, s=rand_ising2d(), β=β_ising2d, niters=10^3, rng=default_rng())
    m, n = size(s)
    prob = ((exp(-2*β*k) for k in -4:4)...,)
    @inbounds for iter in 1:niters, j in 1:n, i in 1:m
        NN = s[ifelse(i == 1, m, i-1), j]
        SS = s[ifelse(i == m, 1, i+1), j]
        WW = s[i, ifelse(j == 1, n, j-1)]
        EE = s[i, ifelse(j == n, 1, j+1)]
        CT = s[i, j]
        k = CT * (NN + SS + WW + EE)
        s[i,j] = ifelse(rand(rng) < prob[k+5], -CT, CT)
    end
    s
end

"""
    ising2d!(::MultiFor, s=rand_ising2d(), β=β_ising2d, niters=10^3, rng=default_rng())

uses the algorithm using multiple for-loops.

Example:


```julia
s = ones_ising2d()
ising2d!(MultiFor(), s)
plot_ising2d(s)
```
"""
function ising2d!(::MultiFor, s=rand_ising2d(), β=β_ising2d, niters=10^3, rng=default_rng())
    m, n = size(s)
    prob = ((exp(-2*β*k) for k in -4:4)...,)
    @inbounds for iter in 1:niters
        let s₁ = s[1,1], k = s₁ * (s[2,1] + s[m,1] + s[1,2] + s[1,n])
            s[1,1] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
        end
        for i in 2:m-1
            let s₁ = s[i,1], k = s₁ * (s[i+1,1] + s[i-1,1] + s[i,2] + s[i,n])
                s[i,1] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
            end
        end
        let s₁ = s[m,1], k = s₁ * (s[1,1] + s[m-1,1] + s[m,2] + s[m,n])
            s[m,1] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
        end
        
        for j in 2:n-1
            let s₁ = s[1,j], k = s₁ * (s[2,j] + s[m,j] + s[1,j+1] + s[1,j-1]) 
                s[1,j] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
            end
            for i in 2:m-1
                let s₁ = s[i,j], k = s₁ * (s[i+1,j] + s[i-1,j] + s[i,j+1] + s[i,j-1])
                    s[i,j] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
                end
            end
            let s₁ = s[m,j], k = s₁ * (s[1,j] + s[m-1,j] + s[m,j+1] + s[m,j-1]) 
                s[m,j] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
            end
        end
        
        let s₁ = s[1,n], k = s₁ * (s[2,n] + s[m,n] + s[1,1] + s[1,n-1])
            s[1,n] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
        end
        for i in 2:m-1
            let s₁ = s[i,n], k = s₁ * (s[i+1,n] + s[i-1,n] + s[i,1] + s[i,n-1])
                s[i,n] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
            end
        end
        let s₁ = s[m,n], k = s₁ * (s[1,n] + s[m-1,n] + s[m,1] + s[m,n-1])
            s[m,n] = ifelse(rand(rng) < prob[k+5], -s₁, s₁)
        end
    end
    s
end

"""
    plot_ising2d(s; size=(201.5, 201.5), color=:gist_earth, clim=(-2, 1.1), kwargs...)

plots the state `s` of 2D Ising model.
"""
function plot_ising2d(s; size=(201.5, 201.5), color=:gist_earth, clim=(-2, 1.1), kwargs...)
    heatmap(s; size=size, color=color, clim=clim,
        colorbar = false, axis = false,
        leftmargin = -20mm, rightmargin = -20mm, 
        top_margin = -2mm, bottom_margin = -20mm,
        titlefontsize = 8 
    )
    plot!(; kwargs...)
end

"""
    gif_ising2d(; s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),
        algorithm=default_algorithm(), progress=true,
        nwarmups=0, nskips=10, niters=100, 
        gifname="ising2d.gif", fps=10,
        size=(201.5, 217.5), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )

creates the gif animation of 2D Ising model.

Example: To create a 200x200 GIF animation, run

```julia
gif_ising2d(s=rand_ising2d(200))
```
"""
function gif_ising2d(; s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),
        algorithm=default_algorithm(), progress=true,
        nwarmups=0, nskips=10, niters=100, 
        gifname="ising2d.gif", fps=10,
        size=(201.5, 217.5), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )
    ising2d!(algorithm, s, β, nwarmups, rng)
    progress && (prog = Progress(niters, 0))
    anim = @animate for t in 0:niters
        iszero(t) || ising2d!(algorithm, s, β, nskips, rng)
        title="β=$(round(β/β_ising2d, digits=4))β_c,  t=$t"
        plot_ising2d(s; size=size, color=color, clim=clim, title=title, kwargs...)
        progress && next!(prog)
    end
    gif(anim, gifname; fps=fps)
end

"""
    mcmc_ising2d!(; s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),
        algorithm=default_algorithm(), progress=true, 
        nwarmups=1000, nskips=100, niters=5000
    )

returns the result of the Markov Chain Monte Carlo simulation with length niters, which is the array of the states of 2D Ising model.
"""
function mcmc_ising2d!(; s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),
        algorithm=default_algorithm(), progress=true, 
        nwarmups=1000, nskips=100, niters=5000
    )
    S = Array{typeof(s), 1}(undef, niters)
    ising2d!(algorithm, s, β, nwarmups, rng)
    progress && (prog = Progress(niters, 1, "mcmc_ising2d!: "))
    for i in 1:niters
        ising2d!(algorithm, s, β, nskips, rng)
        S[i] = copy(s)
        progress && next!(prog)
    end
    S
end

"""
    energy_ising2d(s)

returns the energy density (energy per site) of the state `s` of 2D Ising model.
"""
function energy_density_ising2d(s)
    m, n = size(s)
    E = 0.0
    @inbounds begin
        for j in 1:n, i in 1:m-1
            E -= s[i,j]*s[i+1,j]
        end
        for j in 1:n
            E -= s[m,j]*s[1,j]
        end
        for j in 1:n-1, i in 1:m
            E -= s[i,j]*s[i,j+1]
        end
        for i in 1:m
            E -= s[i,m]*s[i,1]
        end
    end
    E/(m*n)
end

"""
    magnetization_ising2d(s)

returns the magnetization of the state `s` of 2D Ising model.
"""
function magnetization_ising2d(s)
    m, n = size(s)
    sum(s)/(m*n)
end

"""
    plot_mcmc_ising2d(S, E=nothing, M=nothing;
        k=1.0, β=k*β_ising2d, niters=length(S), t=niters,   
        ylim_E=(-1.50, -1.35), ylim_M=(-0.75, 0.75), lw=0.5, alpha=0.8,
        size=(600, 300), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )

plots the MCMC result S with energy per site and magnetization.
"""
function plot_mcmc_ising2d(S, E=nothing, M=nothing;
        k=1.0, β=k*β_ising2d, niters=length(S), t=niters, 
        ylim_E=(-1.55, -1.3), ylim_M=(-0.8, 0.8), lw=0.5, alpha=0.8,
        size=(600, 300), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )
    if isnothing(E)
        E = energy_density_ising2d.(S)
    end
    if isnothing(M)
        M = magnetization_ising2d.(S)
    end
    title="β=$(round(β/β_ising2d, digits=4))β_c,  t=$t"
    P1 = heatmap(S[t]; color=color, clim=clim, colorbar=false, axis=false, 
        title=title)
    P2 = plot(@view(E[1:t]); xlim=(1, niters), ylim=ylim_E, lw=lw, alpha=alpha,
        title="energy per site")
    P3 = plot(@view(M[1:t]); xlim=(1, niters), ylim=ylim_M, lw=lw, alpha=alpha,
        title="magnetization")
    plot(P1, P2, P3; size=size, legend=false, titlefontsize=10,
        layout=@layout([a [b; c]]))
    plot!(; kwargs...)
end

"""
    histogram_mcmc_ising2d(S, E=nothing, M=nothing;
        niters=length(S), bin=min(100, max(10, round(Int, 1.2*√niters))), 
        size=(600, 200), kwargs...
    )

plots the histogram of energy per site and magnetization of the MCMC result S.
"""
function histogram_mcmc_ising2d(S, E=nothing, M=nothing;
        niters=length(S), bin=min(100, max(10, round(Int, 1.2*√niters))), 
        size=(600, 200), kwargs...
    )
    if isnothing(E)
        E = energy_density_ising2d.(S)
    end
    if isnothing(M)
        M = magnetization_ising2d.(S)
    end
    P1 = histogram(E; norm=true, bin=bin, alpha=0.3, title="energy per site")
    P2 = histogram(M; norm=true, bin=bin, alpha=0.3, title="magnetization")
    plot(P1, P2; size=size, legend=false, titlefontsize=10)
    plot!(; kwargs...)
end

"""
    gif_mcmc_ising2d(S=nothing, E=nothing, M=nothing;
        s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),
        algorithm=default_algorithm(), progress=true,
        nwarmups=1000, nskips=1000, niters=500,
        ylim_E=(-1.55, -1.30), ylim_M=(-0.8, 0.8), lw=1.0, alpha=0.8,
        gifname="ising2d_mcmc.gif", fps=10,
        size=(600, 300), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )

creates the gif animation of 2D Ising model with energy per site and magnetization.

Example:
```
gif_mcmc_ising2d()
```
"""
function gif_mcmc_ising2d(S=nothing, E=nothing, M=nothing;
        s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),
        algorithm=default_algorithm(), progress=true,
        nwarmups=1000, nskips=1000, niters=500,
        ylim_E=(-1.55, -1.30), ylim_M=(-0.8, 0.8), lw=1.0, alpha=0.8,
        gifname="ising2d_mcmc.gif", fps=10,
        size=(600, 300), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )
    if isnothing(S)
        S = mcmc_ising2d!(; s=s, k=k, β=β, rng=rng, algorithm=algorithm, 
            progress=progress, nwarmups=nwarmups, nskips=nskips, niters=niters)
    else
        niters = length(S)
    end
    if isnothing(E)
        E = energy_density_ising2d.(S)
    end
    if isnothing(M)
        M = magnetization_ising2d.(S)
    end
    progress && (prog = Progress(niters, 1, "@animate for:  "))
    anim = @animate for t in 1:niters
        plot_mcmc_ising2d(S, E, M; k=k, β=β, t=t,   
            ylim_E=ylim_E, ylim_M=ylim_M, lw=lw, alpha=alpha,
            size=size, color=color, clim=clim, kwargs...
        )
        progress && next!(prog)
    end
    gif(anim, gifname; fps=fps)
end

end # module
