module Ising2D

using Plots
using Plots.PlotMeasures
using Random

export β_ising2d, ising2d!, rand_ising2d, plot_ising2d, gif_ising2d

"""
    β_ising2d = log(1 + √2)/2

is the critical inverse temperature of the 2-dimensional Ising model on the square lattice.
"""
const β_ising2d = log(1 + √2)/2

"""
    ising2d!(s, β=β_ising2d, niters=10^3, rng=Random.default_rng())

updates the 2-dimensional array `s` with values ±1 randomly by the 2-dimensional Ising model rule.
"""
function ising2d!(s, β=β_ising2d, niters=10^3, rng=Random.default_rng())
    m, n = size(s)
    prob = [exp(-2*β*k) for k in -4:4]
    @inbounds for iter in 1:niters, j in 1:n, i in 1:m
        NN = s[ifelse(i == 1, m, i-1), j]
        SS = s[ifelse(i == m, 1, i+1), j]
        WW = s[i, ifelse(j == 1, n, j-1)]
        EE = s[i, ifelse(j == n, 1, j+1)]
        CT = s[i, j]
        k = CT * (NN + SS + WW + EE)
        s[i,j] = ifelse(rand(rng) < prob[k+5], -CT, CT)
    end
end

"""
    rand_ising2d(rng::AbstractRNG, m=100, n=m)
    rand_ising2d(m=100, n=m)

generate the random state of the 2-dimensional Ising model.
"""
rand_ising2d(rng::AbstractRNG, m=100, n=m) = rand(rng, Int8[-1, 1], m, n)
rand_ising2d(m=100, n=m) = rand(Int8[-1, 1], m, n)

macro def(name, definition)
  quote
      macro $(esc(name))()
          esc($(Expr(:quote, definition)))
      end
  end
end

"""
    plot_ising2d(s; 
        size=(201.5, 201.5), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )

plots the state `s` of the 2-dimensional Ising model.
"""
function plot_ising2d(s; 
        size=(201.5, 201.5), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )
    heatmap(s; size=size, color=color, clim=clim,
        colorbar = false, axis = false,
        leftmargin = -20mm, rightmargin = -20mm, 
        top_margin = -2mm, bottom_margin = -20mm,
        titlefontsize = 8 
    )
    plot!(; kwargs...)
end

"""
    gif_ising2d(s = rand_ising2d(), k = 1.0; 
        β = k*β_ising2d, rng = Random.default_rng(),
        nwarmups = 0, nskips = 10, nframes = 100, 
        gifname = "ising2d.gif", fps = 10,
        size=(201.5, 217.5), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )

creates the gif animation of the 2-dimensional Ising model with the initial state `s`.
"""
function gif_ising2d(s = rand_ising2d(), k = 1.0; 
        β = k*β_ising2d, rng = Random.default_rng(),
        nwarmups = 0, nskips = 10, nframes = 100, 
        gifname = "ising2d.gif", fps = 10,
        size=(201.5, 217.5), color=:gist_earth, clim=(-2, 1.1), kwargs...
    )
    ising2d!(s, β, nwarmups, rng)
    anim = @animate for t in 0:nframes
        iszero(t) || ising2d!(s, β, nskips, rng)
        title="β=$(round(β/β_ising2d, digits=4))β_c,  t=$t"
        P = plot_ising2d(s; size=size, color=color, clim=clim, title=title, kwargs...)
    end
    gif(anim, gifname; fps=fps)
end

end # module
