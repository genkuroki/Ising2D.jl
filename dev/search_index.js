var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Ising2D","category":"page"},{"location":"#Ising2D","page":"Home","title":"Ising2D","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Ising2D]","category":"page"},{"location":"#Ising2D.Ising2D","page":"Home","title":"Ising2D.Ising2D","text":"Ising2D is a simple Julia module of 2D Ising model.\n\nBenchmark:\n\nusing Ising2D, BenchmarkTools\ns = rand_ising2d()\nprint(\"VERSION = $VERSION:\")\n@btime ising2d!($s);\n\n\n\n\n\n","category":"module"},{"location":"#Ising2D.β_ising2d","page":"Home","title":"Ising2D.β_ising2d","text":"β_ising2d = log(1 + √2)/2\n\nis the critical inverse temperature of the 2-dimensional Ising model on the square lattice.\n\n\n\n\n\n","category":"constant"},{"location":"#Ising2D.energy_density_ising2d-Tuple{Any}","page":"Home","title":"Ising2D.energy_density_ising2d","text":"energy_ising2d(s)\n\nreturns the energy density (energy per site) of the state s of 2D Ising model.\n\n\n\n\n\n","category":"method"},{"location":"#Ising2D.gif_ising2d-Tuple{}","page":"Home","title":"Ising2D.gif_ising2d","text":"gif_ising2d(; s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),\n    nwarmups=0, nskips=10, niters=100, \n    gifname=\"ising2d.gif\", fps=10,\n    size=(201.5, 217.5), color=:gist_earth, clim=(-2, 1.1), kwargs...\n)\n\ncreates the gif animation of 2D Ising model.\n\nExample: To create a 200x200 GIF animation, run\n\ngif_ising2d(s=rand_ising2d(200))\n\n\n\n\n\n","category":"method"},{"location":"#Ising2D.gif_mcmc_ising2d","page":"Home","title":"Ising2D.gif_mcmc_ising2d","text":"gif_mcmc_ising2d(S=nothing, E=nothing, M=nothing;\n    s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),\n    nwarmups=1000, nskips=1000, niters=500,\n    ylim_E=(-1.55, -1.30), ylim_M=(-0.8, 0.8), lw=1.0, alpha=0.8,\n    gifname=\"ising2d_mcmc.gif\", fps=10,\n    size=(600, 300), color=:gist_earth, clim=(-2, 1.1), kwargs...\n)\n\ncreates the gif animation of 2D Ising model with energy per site and magnetization.\n\nExample:\n\ngif_mcmc_ising2d()\n\n\n\n\n\n","category":"function"},{"location":"#Ising2D.histogram_mcmc_ising2d","page":"Home","title":"Ising2D.histogram_mcmc_ising2d","text":"histogram_mcmc_ising2d(S, E=nothing, M=nothing;\n    niters=length(S), bin=min(100, max(10, round(Int, 1.2*√niters))), \n    size=(600, 200), kwargs...\n)\n\nplots the histogram of energy per site and magnetization of the MCMC result S.\n\n\n\n\n\n","category":"function"},{"location":"#Ising2D.ising2d!","page":"Home","title":"Ising2D.ising2d!","text":"ising2d!(s=rand_ising2d(), β=β_ising2d, niters=10^3, rng=default_rng())\n\nupdates the 2-dimensional array s with values ±1 randomly by the 2-dimensional Ising model rule.\n\nExample: \n\nusing Ising2D\n@time s = ising2d!(rand_ising2d(), β_ising2d, 10^5)\nplot_ising2d(s)\n\n\n\n\n\n","category":"function"},{"location":"#Ising2D.magnetization_ising2d-Tuple{Any}","page":"Home","title":"Ising2D.magnetization_ising2d","text":"magnetization_ising2d(s)\n\nreturns the magnetization of the state s of 2D Ising model.\n\n\n\n\n\n","category":"method"},{"location":"#Ising2D.mcmc_ising2d!-Tuple{}","page":"Home","title":"Ising2D.mcmc_ising2d!","text":"mcmc_ising2d!(; s=rand_ising2d(), k=1.0, β=k*β_ising2d, rng=default_rng(),\n    nwarmups=1000, nskips=100, niters=5000\n)\n\nreturns the result of the Markov Chain Monte Carlo simulation with length niters, which is the array of the states of 2D Ising model.\n\n\n\n\n\n","category":"method"},{"location":"#Ising2D.plot_ising2d-Tuple{Any}","page":"Home","title":"Ising2D.plot_ising2d","text":"plot_ising2d(s; size=(201.5, 201.5), color=:gist_earth, clim=(-2, 1.1), kwargs...)\n\nplots the state s of the 2-dimensional Ising model.\n\n\n\n\n\n","category":"method"},{"location":"#Ising2D.plot_mcmc_ising2d","page":"Home","title":"Ising2D.plot_mcmc_ising2d","text":"plot_mcmc_ising2d(S, E=nothing, M=nothing;\n    k=1.0, β=k*β_ising2d, niters=length(S), t=niters,   \n    ylim_E=(-1.50, -1.35), ylim_M=(-0.75, 0.75), lw=0.5, alpha=0.8,\n    size=(600, 300), color=:gist_earth, clim=(-2, 1.1), kwargs...\n)\n\nplots the MCMC result S with energy per site and magnetization.\n\n\n\n\n\n","category":"function"},{"location":"#Ising2D.rand_ising2d","page":"Home","title":"Ising2D.rand_ising2d","text":"rand_ising2d(rng::AbstractRNG, m=100, n=m)\nrand_ising2d(m=100, n=m)\n\ngenerate the random state of the 2-dimensional Ising model.\n\n\n\n\n\n","category":"function"}]
}
