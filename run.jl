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

iters = setv(params, "iters", int(1e4), int)
num_feat = setv(params, "num_feat", 2, int)
rseed = setv(params, "rseed", rand(Uint), int)
#seed = setv(params, "seed", 1234, int)
seed = setv(params, "seed", rand(Uint), int)
Ntrn = setv(params, "Ntrn", 20, int)
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

kmax = 1
D = 2

######################################################################
# Class 0
######################################################################
#trumu, trucov, dataa, tst_dataa,lam1,lam2 = gen_data_jason(0.0)
#tst_data = rand(510,D) .* 500 .+ 8
#tst_data[:,1] = 0
##tst_data[:,1] = rand(0:1, 510)
#dataa = tst_data[1:10,:]


######################################################################
# Class 1
######################################################################
#trumub, trucovb, datab, tst_datab,_,_ = gen_data_jason(-1.0)
#tst_datab = rand(510,D) .* 500 .+ 8
#tst_datab[:,1] = 2
##tst_datab[:,1] = rand(0:1, 510)
#datab = tst_datab[1:10,:]

errs = Float64[]
dtries = logspace(0,3,40)
for d=dtries

    cls = MPM.mpm_classifier(dataa, datab; burn=1000, thin=50, d=d, usepriors=false)
    @time MPM.sample(cls, 10000)
    #pts1 = MPM.gen_posterior_points(100, d, cls.mcmc1.db)
    #pts2 = MPM.gen_posterior_points(100, d, cls.mcmc2.db)

    #@time beis,r = MPM.bee_e_nsum(cls, 50)
    @time bemc = MPM.bee_e_mc(cls, dmean=d)
    @show bemc
    push!(errs, bemc)

    #N = 30
    #beis = Float64[]
    #for i=1:N
        #b,r = MPM.bee_e_nsum(cls, 50)
        #push!(beis, b[3])
    #end
    #@show mean(beis), std(beis)
    #beis = [MPM.bee_e_mc(cls, numpts=50) for i=1:N]
    #@show mean(beis), std(beis)

end

#err = MPM.error_points(mymh_a.db, mymh_b.db, [tst_data; tst_datab], [zeros(size(tst_data,1)), ones(size(tst_datab,1))]) 
#println("holdout error: $err")

######################################################################
# Plotting
######################################################################

#include("src/plot_utils.jl")
#mins, maxs, (lens, steps, grid) = MPM.get_grid(vcat(data,datab))

#close("all")

#ga = MPM.calc_g(grid, mymh_a.db, 20)
#ga = exp(ga.-maximum(ga))
#imshow(reshape(ga,lens...)', extent=[mins[1],maxs[1],mins[2],maxs[2]], aspect="equal", origin="lower")
#colorbar()
#contour(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
#plot(data[:,1], data[:,2], "g.", alpha=0.8)

#allnodes = collect(r)
#for n in allnodes
    #if length(n.vals)!=0
        #plot(n.mins[1], n.mins[2], "g.")
    #end
#end

#figure()
#gavg = ga - gb
#imshow(reshape(gavg,n1,n2)', extent=gext, aspect="equal", origin="lower")
#colorbar()
#contour(reshape(gavg,n1,n2)', (0.0,), extent=gext, aspect="equal", origin="lower")
#plot(data[:,1], data[:,2], "g.", alpha=0.8)
#plot(datab[:,1], datab[:,2], "r.", alpha=0.8)

reload(Pkg.dir("OBC","src","plot_utils.jl"))
fa() = plot_traces(cls.mcmc1.db, [:mu,:sigma,:lam,:energy])
