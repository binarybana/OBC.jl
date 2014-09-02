reload("OBC.jl")
#using PyPlot
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

iters = setv(params, "iters", int(3e3), int)
num_feat = setv(params, "num_feat", 2, int)
rseed = setv(params, "rseed", rand(Uint), int)
#seed = setv(params, "seed", 1234, int)
seed = setv(params, "seed", rand(Uint), int)
Ntrn = setv(params, "Ntrn", 10, int)
Ntst = setv(params, "Ntst", 500, int)
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
    ps = zeros(Int, size(lams))
    for i in 1:size(lams,1)
        for j in 1:size(lams,2)
            ps[i,j] = rand(Poisson(rand(Uniform(lowd,highd)) * exp(lams[i,j])))
        end
    end
    lams, ps
end

function gen_data_jason(mu)
    D = num_feat
    lmu = zeros(D) .+ mu 
    cov = eye(D)*0.5
    cov[1,2] = -0.1
    cov[2,1] = -0.1
    lam1, trn_data = gen_data(lmu, cov, Ntrn)
    lam2, tst_data = gen_data(lmu, cov, Ntst)
    return lmu, cov, trn_data', tst_data', lam1, lam2
end

D = num_feat

######################################################################
# Class 0
######################################################################
trumu, trucov, dataa, tst_dataa,lam1,lam2 = gen_data_jason(0.0)
#tst_data = rand(510,D) .* 500 .+ 8
#tst_data[:,1] = 0
##tst_data[:,1] = rand(0:1, 510)
#dataa = tst_data[1:10,:]


######################################################################
# Class 1
######################################################################
trumub, trucovb, datab, tst_datab,_,_ = gen_data_jason(-1.0)
#tst_datab = rand(510,D) .* 500 .+ 8
#tst_datab[:,1] = 2
##tst_datab[:,1] = rand(0:1, 510)
#datab = tst_datab[1:10,:]

#d=100.0
d1 = ones(Ntrn)*100.0 .+ randn(Ntrn)*10
d2 = ones(Ntrn)*100.0 .+ randn(Ntrn)*10
dataa .*= d1 ./ 100.0
datab .*= d2 ./ 100.0
dmean1,dmean2 = mean(d1),mean(d2)
dmean = mean([dmean1,dmean2])

cls = MPM.mpm_classifier(dataa, datab, burn=1000, thin=50, d1=d1, d2=d2, kappa=50.0, usepriors=true)
@time MPM.sample(cls, 10000)#, verbose=true)
@show bemc = MPM.bee_e_mc(cls, (dmean1,dmean2),numpts=50)
@show bemc = MPM.bee_e_mc(cls, (dmean1,dmean2),numpts=50)
@show beem1,beem2 = MPM.bee_moments(cls, (dmean1,dmean2))
@show beem1,beem2 = MPM.bee_moments(cls, (dmean1,dmean2))
@show MPM.bee_moments_2sample(cls, (dmean1,dmean2),numpts=50)
@show MPM.bee_moments_2sample(cls, (dmean1,dmean2),numpts=50)

#function project!(obj::MPM.MPMParams, sel)
    #obj.mu = obj.mu[sel]
    #obj.sigma = obj.sigma[sel,sel]
    #obj.sigpre = obj.sigpre[sel,sel]
    #obj.lam = obj.lam[sel,:]
    #return obj
#end

#function project(cls::OBC.BinaryClassifier,sel)
    #subcls = deepcopy(cls)
    #map!(x->project!(x,sel), subcls.mcmc1.db)
    #map!(x->project!(x,sel), subcls.mcmc2.db)
    #project!(subcls.cls1.curr,sel)
    #project!(subcls.cls1.old,sel)
    #project!(subcls.cls2.curr,sel)
    #project!(subcls.cls2.old,sel)
    #return subcls
#end

#for sel = combinations(1:3,2)
    #subcls = project(cls,sel)
    #@show sel
    #@show bemc = MPM.bee_e_mc(subcls, (dmean1,dmean2))
#end

#err = MPM.error_points(cls, [tst_dataa; tst_datab], [zeros(size(tst_dataa,1)), ones(size(tst_datab,1))], dmean=dmean) 
#println("holdout error: $err")

#cls = MPM.mpm_classifier(dataa, datab, burn=1000, thin=50, d1=dmean, d2=dmean, usepriors=false)
#@time MPM.sample(cls, 10000)#, verbose=true)
#@show bemc = MPM.bee_e_mc(cls, dmean=dmean)
#@show beem1,beem2 = MPM.bee_moments(cls, dmean=dmean)

#for i=1:2
    #@time bemc = MPM.bee_e_mc(cls, dmean=dmean)
    #@time beem1,beem2 = MPM.bee_moments(cls, dmean=dmean)

    #@show bemc
    #@show beem1, beem2
    #@show beemse = beem2 - beem1^2
#end

#err = MPM.error_points(cls, [tst_dataa; tst_datab], [zeros(size(tst_dataa,1)), ones(size(tst_datab,1))], dmean=dmean) 
#println("holdout error: $err")

