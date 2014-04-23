reload("src/OBC.jl")
reload("src/mpm.jl")
#using OBC, MPM
using PyPlot
using Distributions

#srand(1)
blas_set_num_threads(1)


#if 'PARAM' in os.environ:
    #params = yaml.load(os.environ['PARAM'])
#else:
    #params = {}
params = (ASCIIString=>Any)[]

function setv(p,s,d,conv=None)
    if !(s in p)
    #if !contains(p,s)
        p[s] = d # Possible error TODO str(d)
        return d
    elseif conv != None
        return conv(p[s])
    else
        return p[s]
    end
end

iters = setv(params, "iters", int(1e4), int)
num_feat = setv(params, "num_feat", 2, int)
rseed = setv(params, "rseed", rand(Uint), int)
#seed = setv(params, "seed", 1234, int)
seed = setv(params, "seed", rand(Uint), int)
Ntrn = setv(params, "Ntrn", 10, int)
Ntst = setv(params, "Ntst", 300, int)
f_glob = setv(params, "f_glob", 1, int)
subclasses = setv(params, "subclasses", 2, int)
f_het = setv(params, "f_het", 1, int)
f_rand = setv(params, "f_rand", 0, int)
rho = setv(params, "rho", 0.6, float)
f_tot = setv(params, "f_tot", f_glob+f_het*subclasses+f_rand, int)
blocksize = setv(params, "blocksize", 1, int)
mu0 = setv(params, "mu0", 0.2, float)
mu1 = setv(params, "mu1", 0.6, float)
sigma0 = setv(params, "sigma0", 0.1, float)
sigma1 = setv(params, "sigma1", 0.5, float)
lowd = setv(params, "lowd", 9.0, float)
highd = setv(params, "highd", 11.0, float)
numlam = setv(params, "numlam", 20, int)

function gen_data(mu, cov, n)
    lams = rand(MultivariateNormal(mu, cov), n)
    #println("lam stats:")
    #println("$(maximum(lams)) $(minimum(lams)) $(mean(lams)) $(std(lams))")
    ps = zeros(size(lams))
    for i in 1:size(lams,1)
        for j in 1:size(lams,2)
            ps[i,j] = rand(Poisson(rand(Uniform(lowd,highd)) * exp(lams[i,j])))
        end
    end
    ps
end

function gen_data_jason(mu)
    D = num_feat
    lmu0 = zeros(D) .+ mu 
    cov0 = eye(D)*0.05
    #@show lmu0 cov0
    trn_data0 = gen_data(lmu0, cov0, Ntrn)
    tst_data0 = gen_data(lmu0, cov0, Ntst)
    #return lmu0, cov0, ones(2,10)+1, tst_data0
    return lmu0, cov0, trn_data0', tst_data0'
end

kmax = 1
D = 2

######################################################################
# Class 0
######################################################################

#trumu, trucov, data, tst_data = gen_data_jason(0.0)
tst_data = rand(510,D) .+ 10
#tst_data[:,1] .+= 0
tst_data[:,1] = rand(0:1, 510)
data = tst_data[1:10,:]

cov = eye(D,D) 
#cov = [0.5 -0.2; -0.2 1.0]
mu = ones(D,kmax)
#mu = [1.6 -1.0; 1.6 -1.0]
#mu = [0.2 0.2]'
prior = MPM.MPMPrior(D=D, kmax=kmax, kappa=10.)
start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
    cov, #sigma :: Matrix{Float64}
    ones(kmax)/kmax, #w :: Vector{Float64}
    clamp(log(data'/10),-8.0,Inf), #lam :: Matrix{Float64}
    1) #k :: Int

pmoves = MPM.MPMPropMoves()
pmoves.lammove = 0.01
pmoves.mumove = 0.1
pmoves.priorkappa = 100.0

obj_a = MPM.MPMCls(prior, data, deepcopy(start), pmoves, 10.0)
obj_a.usepriors = true

mymh_a = OBC.MHRecord(obj_a,burn=1000,thin=50)
sample(mymh_a,10000)

######################################################################
# Class 1
######################################################################

#trumub, trucovb, datab, tst_datab = gen_data_jason(-1.0)
tst_datab = rand(510,D) .+ 10
tst_datab[:,2] .+= 0.0
tst_datab[:,1] = rand(0:1, 510)
datab = tst_datab[1:10,:]

start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
    cov, #sigma :: Matrix{Float64}
    ones(kmax)/kmax, #w :: Vector{Float64}
    clamp(log(datab'/10),-8.0,Inf), #lam :: Matrix{Float64}
    1) #k :: Int
obj_b = MPM.MPMCls(prior, datab, deepcopy(start), pmoves, 10.0)

mymh_b = OBC.MHRecord(obj_b,burn=1000,thin=50)
OBC.sample(mymh_b,10000)

#gext = [0,30.0,0,30]
#n1, n2, xstride, ystride, grid = MPM.gen_grid(gext, 30)
#gext, (n1, n2, xstride, ystride, grid) = MPM.get_grid(vcat(data,datab))

#g1 = MPM.calc_g(grid, mymh_a.db, 20)
#g2 = MPM.calc_g(grid, mymh_b.db, 20)

@time be1 = MPM.bee_e_data(vcat(data,datab), mymh_a.db, mymh_b.db, 20)
@time be2 = MPM.bee_e_data(vcat(data,datab), mymh_a.db, mymh_b.db, 20)
@time becube = MPM.bee_e_cube(vcat(data,datab), mymh_a.db, mymh_b.db, 20)
@show be1
@show be2
@show becube

err = MPM.error_points(mymh_a.db, mymh_b.db, [tst_data; tst_datab], [zeros(size(tst_data,1)), ones(size(tst_datab,1))]) 
println("holdout error: $err")

######################################################################
# Plotting
######################################################################

#include("src/plot_utils.jl")
#n1, n2, gext, grid = get_grid(tst_data)
#n1, n2, gext, grid = gen_grid([0,10.0,0,10])

#gavg,lams = calc_g(grid, [start], 1000)
#imshow(reshape(gavg,N,N), extent=gext, aspect="equal", origin="lower")
#colorbar()
#plot(data[:,1], data[:,2], "g.", alpha=0.8)
#plot(10*exp(lams[1,:]), 10*exp(lams[2,:]), "k.")

#close("all")

#ga = calc_g(grid, mymh_a.db, 20)
#imshow(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
#colorbar()
#contour(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
#plot(data[:,1], data[:,2], "g.", alpha=0.8)

#figure()
#gb = calc_g(grid, mymh_b.db, 20)
#imshow(reshape(gb,n1,n2)', extent=gext, aspect="equal", origin="lower")
#colorbar()
#contour(reshape(gb,n1,n2)', extent=gext, aspect="equal", origin="lower")
#plot(datab[:,1], datab[:,2], "m.", alpha=0.8)

#figure()
#gavg = ga - gb
#imshow(reshape(gavg,n1,n2)', extent=gext, aspect="equal", origin="lower")
#colorbar()
#contour(reshape(gavg,n1,n2)', (0.0,), extent=gext, aspect="equal", origin="lower")
#plot(data[:,1], data[:,2], "g.", alpha=0.8)
#plot(datab[:,1], datab[:,2], "r.", alpha=0.8)

#gmin = min(gg0,gg1)
#close("all")
#figure()
#N=n1
#imshow(reshape(gmin,N,N)', extent=gext, aspect="equal", origin="lower")
#contour(reshape(gmin,N,N)', extent=gext, aspect="equal", origin="lower")
#plot(data[:,1], data[:,2], "g.", alpha=0.8)
#plot(datab[:,1], datab[:,2], "r.", alpha=0.8)

#figure()
#gmin = exp(min(g0,g1))
#imshow(reshape(gmin,N,N)', extent=gext, aspect="equal", origin="lower")
#contour(reshape(gmin,N,N)', extent=gext, aspect="equal", origin="lower")
#plot(data[:,1], data[:,2], "g.", alpha=0.8)
#plot(datab[:,1], datab[:,2], "r.", alpha=0.8)

#figure()
#gmin = g0 .- g1
#imshow(reshape(gmin,N,N)', extent=gext, aspect="equal", origin="lower")
#contour(reshape(gmin,N,N)', levels = [0.0], extent=gext, aspect="equal", origin="lower")
#plot(tst_data[:,1], tst_data[:,2], "g.", alpha=0.8)
#plot(tst_datab[:,1], tst_datab[:,2], "r.", alpha=0.8)

#@time g0o = MPM.calc_g(grid, mymh_a.db, 20)
#@time g1o = MPM.calc_g(grid, mymh_b.db, 20)
#@show be2, be2test = MPM.e_error_eff(g0o, g1o, xstride*ystride)

#g() = plot(tst_data[:,1]+rand(size(tst_data,1)), tst_data[:,2]+rand(size(tst_data,1)), "g.", alpha=0.3)
#fa() = plot_traces(mymh_a.db, [:mu,:sigma,:lam,:energy])
#fb() = plot_traces(mymh_b.db, [:mu,:sigma,:lam,:energy])
