# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.2'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Julia 1.6.0-beta1 depwarn -O3
#     language: julia
#     name: julia-1.6-depwarn-o3
# ---

# %% [markdown]
# # Ising2D
#
# * Copyright (c) 2020 Gen Kuroki
# * License: MIT
# * Repository: https://github.com/genkuroki/Ising2D.jl

# %% [markdown] {"toc": true}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Benchmarks" data-toc-modified-id="Benchmarks-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Benchmarks</a></span></li><li><span><a href="#Examples" data-toc-modified-id="Examples-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Examples</a></span></li><li><span><a href="#Documents" data-toc-modified-id="Documents-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Documents</a></span></li></ul></div>

# %%
if isfile("Project.toml")
    using Pkg
    Pkg.activate(".")
    using Revise
end

# %%
using Ising2D, Plots, Random

# %% [markdown]
# ## Benchmarks

# %%
using Ising2D, Random
Random.seed!(464937337564)
print("IfElse() on v$VERSION:")
@time s = ising2d!(rand_ising2d(), β_ising2d, 10^5; algorithm=IfElse())
plot_ising2d(s)

# %%
using Ising2D, Random
Random.seed!(464937337564)
print("MultiFor() on v$VERSION:")
@time s = ising2d!(rand_ising2d(), β_ising2d, 10^5; algorithm=MultiFor())
plot_ising2d(s)

# %%
using Ising2D, BenchmarkTools, Random
Random.seed!(464937337564)
s = rand_ising2d()
print("IfElse() on v$VERSION:")
@btime ising2d!($s; algorithm=IfElse());

# %%
using Ising2D, BenchmarkTools, Random
Random.seed!(464937337564)
s = rand_ising2d()
print("MultiFor() on v$VERSION:")
@btime ising2d!($s; algorithm=MultiFor());

# %% [markdown]
# ## Examples

# %%
Random.seed!(464937337564);

# %%
s = rand_ising2d(200)
P0 = plot_ising2d(s)

# %%
ising2d!(s, β_ising2d, 500)
P1 = plot_ising2d(s)

# %%
png(P0, "images/s0.png")
png(P1, "images/s1.png")

# %%
s = rand_ising2d(200)
gif_ising2d(s=s, k=1.0, nwarmups=0, nskips=1, niters=500, fps=15,
    gifname="images/ising2d.gif")

# %%
s₀ = rand_ising2d()
@show energy_density_ising2d(s₀)
@show magnetization_ising2d(s₀)
s = copy(s₀)
S = mcmc_ising2d!(s=s, k=1.0, nwarmups=1000, nskips=100, niters=5000)
E = energy_density_ising2d.(S)
M = magnetization_ising2d.(S)
@show mean(E)
@show mean(M)

histogram_mcmc_ising2d(S)

# %%
plot_mcmc_ising2d(S)

# %%
plot_mcmc_ising2d(S; t=4000)

# %%
gif_mcmc_ising2d(S[10:10:end]; gifname="images/ising2d_mcmc.gif")

# %% [markdown]
# ## Documents

# %%
?Ising2D

# %%
?Ising2D.β_ising2d

# %%
?Ising2D.rand_ising2d

# %%
?Ising2D.ones_ising2d

# %%
?Ising2D.ising2d!

# %%
?Ising2D.plot_ising2d

# %%
?Ising2D.gif_ising2d

# %%
?Ising2D.mcmc_ising2d!

# %%
?Ising2D.energy_density_ising2d

# %%
?Ising2D.magnetization_ising2d

# %%
?Ising2D.plot_mcmc_ising2d

# %%
?Ising2D.histogram_mcmc_ising2d

# %%
?Ising2D.gif_mcmc_ising2d

# %%
