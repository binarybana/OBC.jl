using PyPlot
using Distributions

close("all")
include("plot_utils.jl")
N, gext, grid = gen_grid([0,1000.0,0,1000])

x = sum(logpdf(Poisson(10*exp(4)), grid), 2)

imshow(reshape(x,N,N)', extent=gext, aspect="equal", origin="lower")
colorbar()
contour(reshape(x,N,N)', extent=gext, aspect="equal", origin="lower")

