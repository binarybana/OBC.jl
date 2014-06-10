reload("OBC")
using HDF5, JLD
using PyPlot
using Distributions

kmax = 1
D = 2
iters = 10000

loc = "/home/bana/largeresearch/rnaseq.jld"
if !isfile(loc)
    loc = "/home/bana/largeresearch/seq-data/rnaseq.jld"
end
@load loc brca_norm brca_norm_sampids brca_norm_geneids
@load loc luad_norm luad_norm_sampids luad_norm_geneids

blas_set_num_threads(1)

#pick two random genes
#data0 = brca_norm[[5,4,3],[1:40]]'
#data1 = luad_norm[[5,4,3],[1:40]]'
data0 = brca_norm[[5,3],[1:40]]'
data1 = luad_norm[[5,3],[1:40]]'
#data0 = convert(Matrix{Int}, round(brca_norm[:,[2,3]]))
#data1 = convert(Matrix{Int}, round(luad_norm[:,[2,3]]))

#println(brca_norm_geneids[[2,3]])
#println(luad_norm_geneids[[2,3]])

#plot(data0[:,1], data0[:,2], "go", label="brca")
#plot(data1[:,1], data1[:,2], "ro", label="luad")

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
prior = MPM.MPMPrior(D=D, kmax=kmax, kappa=10., S=eye(2)*10)
start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
    cov, #sigma :: Matrix{Float64}
    ones(kmax)/kmax, #w :: Vector{Float64}
    clamp(log(data0'/d),-8.0,Inf), #lam :: Matrix{Float64}
    1) #k :: Int

pmoves = MPM.MPMPropMoves()
pmoves.lammove = 0.005
pmoves.mumove = 0.6
pmoves.priorkappa = 100.0

obj_a = MPM.MPMSampler(prior, data0, deepcopy(start), pmoves, d)
obj_a.usepriors = true

mymh_a = MPM.MHRecord(obj_a,burn=1000,thin=50)
sample(mymh_a,iters)

######################################################################
# Class 1
######################################################################

start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
    cov, #sigma :: Matrix{Float64}
    ones(kmax)/kmax, #w :: Vector{Float64}
    clamp(log(data1'/d),-8.0,Inf), #lam :: Matrix{Float64}
    1) #k :: Int
obj_b = MPM.MPMSampler(prior, data1, deepcopy(start), pmoves, d)
mymh_b = MPM.MHRecord(obj_b,burn=1000,thin=50)
sample(mymh_b,iters)

######################################################################
# BEE Computation
######################################################################

#@time be1 = MPM.bee_e_data(vcat(data0,data1), mymh_a.db, mymh_b.db, 20)
#@time be2 = MPM.bee_e_data(vcat(data1,data1), mymh_a.db, mymh_b.db, 20)
@time bec = MPM.bee_e_cube(vcat(data1,data1), mymh_a.db, mymh_b.db, 20)

#@show be1
#@show be2
@show bec

assert(false)
######################################################################
# Plotting
######################################################################

include("src/plot_utils.jl")
mins, maxs, (lens, steps, grid) = MPM.get_grid(vcat(data0, data1), factor=2, N=100)
#n1, n2, gext, grid = gen_grid([0,10.0,0,10])

#gavg,lams = MPM.calc_g(grid, [start], 1000)
#imshow(reshape(gavg,N,N), extent=gext, aspect="equal", origin="lower")
#colorbar()
#plot(data0[:,1], data0[:,2], "g.", alpha=0.8)
#plot(10*exp(lams[1,:]), 10*exp(lams[2,:]), "k.")

#close("all")

ga = MPM.calc_g(grid, mymh_a.db, 20, dmean=d)
#imshow(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
#colorbar()
#contour(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
#plot(data0[:,1], data0[:,2], "g.", alpha=0.8)

#figure()
gb = MPM.calc_g(grid, mymh_b.db, 20, dmean=d)
#imshow(reshape(gb,n1,n2)', extent=gext, aspect="equal", origin="lower")
#colorbar()
#contour(reshape(gb,n1,n2)', extent=gext, aspect="equal", origin="lower")
#plot(data1[:,1], data1[:,2], "m.", alpha=0.8)

figure()
gavg = exp(min(ga,gb))
#gavg = ga - gb
#gavg = ga
gext = [mins[1], maxs[1], mins[2], maxs[2]]
imshow(reshape(gavg,lens[1],lens[2])', extent=gext, aspect="equal", origin="lower", interpolation="nearest")
colorbar()
#contour(reshape(gavg,lens[1],lens[2])', (0.0,), extent=gext, aspect="equal", origin="lower")
plot(data0[:,1], data0[:,2], "g.", alpha=0.8)
plot(data1[:,1], data1[:,2], "r.", alpha=0.8)

#println("$(MPM.error_points(mymh_a.db, mymh_b.db, [tst_data; tst_data1], [zeros(size(tst_data,1)), ones(size(tst_data1,2))] ))")

#g() = plot(tst_data[:,1]+rand(size(tst_data,1)), tst_data[:,2]+rand(size(tst_data,1)), "g.", alpha=0.3)
#fa() = plot_traces(mymh_a.db, [:mu,:sigma,:lam,:energy])
#fb() = plot_traces(mymh_b.db, [:mu,:sigma,:lam,:energy])
