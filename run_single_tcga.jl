using OBC, MPM
using HDF5, JLD
using PyPlot
using Distributions
using Stats

loc = "/home/bana/largeresearch/rnaseq.jld"
if !isfile(loc)
    loc = "/home/bana/largeresearch/seq-data/rnaseq.jld"
end
@load loc luad_norm luad_norm_sampids luad_norm_geneids

d = 10.0
kmax = 1
iters = 15000
burn = 2000
thin = 50
D = 2
rounds = 1
numsamps = 100

Ts = [x->quantile(x,[0.05])[1], x->quantile(x,[0.5])[1], x->quantile(x,[0.95])[1], iqr, var]
Tdesc = ["5th quantile", "Median", "95th quantile", "IQR", "Variance"]

pvals = zeros(D*rounds, length(Ts))
true_ts = similar(pvals)
tvals = zeros(D*rounds, length(Ts), int((iters-burn)/thin + 1 ))

samps = sample(1:size(luad_norm,2), numsamps)
genes = sample(1:size(luad_norm,1), D*rounds)

######################################################################
# Fit Model
######################################################################
cov = eye(D,D) 
#cov = [0.5 -0.2; -0.2 1.0]
mu = ones(D,kmax)
#mu = [1.6 -1.0; 1.6 -1.0]
#mu = [0.2 0.2]'
prior = MPM.MPMPrior(D=2, kmax=kmax, kappa=10., S=eye(2)*10)

pmoves = MPM.MPMPropMoves()
pmoves.lammove = 0.005
pmoves.mumove = 0.08
pmoves.priorkappa = 100.0
global data, mymh_a

for i=1:rounds

    #pick two random genes
    data = luad_norm[genes[(i-1)*D+1:i*D],samps]'

    start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
        cov, #sigma :: Matrix{Float64}
        ones(kmax)/kmax, #w :: Vector{Float64}
        clamp(log(data'/d),-10.0,Inf), #lam :: Matrix{Float64}
        1) #k :: Int

    obj_a = MPM.MPMSampler(prior, data, deepcopy(start), pmoves, d)
    obj_a.usepriors = true

    mymh_a = MPM.MHRecord(obj_a,burn=burn,thin=thin)
    @time sample(mymh_a,iters)
    @time temp1,temp2,temp3 = MPM.calc_pvals(Ts, data', mymh_a.db, d)
    pvals[(i-1)*D+1:i*D, :] = temp1
    true_ts[(i-1)*D+1:i*D, :] = temp2
    tvals[(i-1)*D+1:i*D, :, :] = temp3
end

#assert(false)
######################################################################
# Plotting
######################################################################

include("src/plot_utils.jl")
n1, n2, gext, grid = get_grid(data)
#n1, n2, gext, grid = gen_grid([0,10.0,0,10])

close("all")

#@time ga,lams = MPM.calc_g(grid, mymh_a.db, 20, dmean=d)
#imshow(reshape(ga,n1,n2)', extent=gext, origin="lower")
#imshow(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
#colorbar()
#contour(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
fig = figure()
ax1 = fig[:add_subplot](211)
plot(data[:,1], data[:,2], "g.", alpha=0.5)
plt.grid(true)

title("TCGA Data")
newpts = MPM.gen_posterior_points(size(data,1), d, mymh_a.db) 
newpts = newpts[:,sample(1:size(newpts,2), size(data,1))]
fig[:add_subplot](212, sharex=ax1, sharey=ax1)
title("Replicated data from predictive posterior")
plt.grid(true)
plot(newpts[1,:], newpts[2,:], "r.", alpha=0.5)

xlabel(string(luad_norm_geneids[genes[1]], " counts"))
#ylabel(string(luad_norm_geneids[genes[2]], " counts"))
fig[:text](0.06, 0.5, string(luad_norm_geneids[genes[2]], " counts"), ha="center", va="center", rotation="vertical")

function plot_hist(geneid, tid)
    figure()
    plt.hist(vec(tvals[geneid,tid,:]), bins=30)
    axvline(true_ts[geneid,tid], color="k", linewidth=4)
end

#function plot_pval(q)
    #p,trut,tvals = MPM.calc_pval(y->quantile(vec(y[1,:]),[q])[1], data', mymh_a.db, d)
    #figure()
    #plt.hist(tvals, bins=30)
    #axvline(trut, color="k", linewidth=2)
    #@show q,p
#end

#plot_pval(0.05)
#plot_pval(0.5)
#plot_pval(0.95)

#g() = plot(tst_data[:,1]+rand(size(tst_data,1)), tst_data[:,2]+rand(size(tst_data,1)), "g.", alpha=0.3)
#fa() = plot_traces(mymh_a.db, [:mu,:sigma,:lam,:energy])

#for i in 1:size(all,1)
       #for j in 1:size(all,2)
       #if j==1
       #print(all[i,j], " & ")
       #else
       #@printf "%.2f %s" all[i,j] " & "
       #end
       #end
       #println(" \\\\")
       #end
