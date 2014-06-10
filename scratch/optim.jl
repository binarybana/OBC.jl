using OBC, MPM
using HDF5, JLD
using PyPlot
using Distributions

kmax = 1
D = 2
iters = 5000

loc = "/home/bana/largeresearch/rnaseq.jld"
if !isfile(loc)
    loc = "/home/bana/largeresearch/seq-data/rnaseq.jld"
end
@load loc luad_norm luad_norm_sampids luad_norm_geneids

#pick two random genes
data = luad_norm[[4,5],:]'
#data = luad_norm[[2,3],:]'

#println(luad_norm_geneids[[2,3]])
#plot(data[:,1], data[:,2], "go", label="brca")
#assert(false)

d = 100.0
######################################################################
# Class 0
######################################################################
cov = eye(D,D) 
#cov = [0.5 -0.2; -0.2 1.0]
mu = ones(D,kmax)
#mu = [1.6 -1.0; 1.6 -1.0]
#mu = [0.2 0.2]'
prior = MPM.MPMPrior(D=2, kmax=kmax, kappa=10., S=eye(2)*10)
start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
    cov, #sigma :: Matrix{Float64}
    ones(kmax)/kmax, #w :: Vector{Float64}
    clamp(log(data'/d),-3.0,Inf), #lam :: Matrix{Float64}
    1) #k :: Int

pmoves = MPM.MPMPropMoves()
obj_a = MPM.MPMSampler(prior, data, deepcopy(start), pmoves, d)
obj_a.usepriors = true

opt = Opt(:LN_COBYLA, 2)
lower_bounds!(opt, [-Inf, 0.])
xtol_rel!(opt,1e-4)

min_objective!(opt, myfunc)
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,2,0), 1e-8)
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,-1,1), 1e-8)

(minf,minx,ret) = optimize(opt, [1.234, 5.678])
println("got $minf at $minx after $count iterations (returned $ret)")
