module MPM
    
using Stats
using Distributions
using OBC
using Cubature
using Iterators
import OBC: propose, energy, reject
include("utils.jl")

export MPMPrior, MPMParams, MPMPropMoves, MPMCls, calc_g, gen_points, 
    gen_posterior_points, calc_pvals, error_points, predict, error_moments_cube, 
    var_error_hists, e_error_hists, e_error_eff

######################################################################
######################################################################
######################################################################
######################################################################

type MPMPrior
    D :: Int
    mu :: Matrix{Float64}
    mu_sigmas :: Vector{Float64} # Unvariate variances on mu
    kappa :: Float64
    S :: Matrix{Float64}
    mix_comps :: Float64
    kmax :: Int
end

function MPMPrior(;D=2, kmax=2, args...) 
    kappa = 6.0
    p = MPMPrior(D, zeros(D,kmax), ones(D)*30, kappa, (kappa-1-D)*eye(D,D), 0.8, kmax)
    for (sym, val) in args
        if sym in MPMPrior.names
            setfield!(p, sym, val)
        else
            error("Quantity $sym not recognized in MPMPrior constructor")
        end
    end
    p
end

######################################################################
######################################################################
######################################################################
######################################################################

type MPMParams 
    mu :: Matrix{Float64}
    sigma :: Matrix{Float64}
    w :: Vector{Float64}
    lam :: Matrix{Float64}
    k :: Int
    energy :: Float64 # I dunno... does this belong here?
    # d :: Vector{Float64}
    function MPMParams(mu, sigma, w, lam, k)
        @assert size(mu,2) == k
        @assert size(sigma) == (size(mu,1),size(mu,1))
        @assert length(w) == k
        @assert size(lam,1) == size(mu,1)
        new(mu, sigma, w, lam, k, 0.0)
    end
end

function copy!(to::MPMParams, from::MPMParams)
    to.mu[:] = from.mu
    to.sigma[:] = from.sigma
    to.w[:] = from.w
    to.lam[:] = from.lam
    to.k = from.k
    to.energy = from.energy
end

#draw_params(p::MPMPrior) = MPMParams(rand(Normal(p.theta_mu, p.theta_sigma)),

#function draw_data(params::MPMParams, n=10)
#end

######################################################################
######################################################################
######################################################################
######################################################################

type MPMPropMoves
    mumove :: Float64
    lammove :: Float64
    wmove :: Float64
    birthmove :: Float64
    priorkappa :: Float64
end
MPMPropMoves() = MPMPropMoves(0.1, 0.01, 0.1, 0.1, 80.0)


type MPMCls <: Sampler
    curr :: MPMParams
    old :: MPMParams
    prior :: MPMPrior
    data :: Matrix{Float64}
    propmoves :: MPMPropMoves
    usepriors :: Bool
    green_factor :: Float64
    d :: Vector{Float64}
end

function MPMCls(prior::MPMPrior, data::Matrix{Float64}, obj::MPMParams, propmoves::MPMPropMoves, d::Float64)
    @assert size(data,1) > size(data,2)
    MPMCls(deepcopy(obj), deepcopy(obj), prior, data, propmoves, 
            true, 0.0, ones(size(data,1))*d) # FIXME hardcoded d
end

function propose(obj::MPMCls)
    curr = obj.curr
    prior = obj.prior
    propmoves = obj.propmoves

    copy!(obj.old, curr)

    range = [0:2]
    obj.green_factor = 0.0
    
    if curr.k == 1 && curr.k == prior.kmax
        scheme = 2
    elseif curr.k == 1
        scheme = wsample(range, [1., 0., 7.])
    elseif curr.k == prior.kmax
        scheme = wsample(range, [0., 1., 7.])
    else
        scheme = wsample(range, [1., 1., 6.])
    end
    scheme = 2

    if scheme == 0 # Birth
        curr.mu[:,curr.k+1] = curr.mu[:,curr.k]
        addfac = randn(prior.D) * propmoves.birthmove
        obj.green_factor -= sum(logpdf(Normal(0.0, propmoves.birthmove), addfac))
        curr.mu[:,curr.k+1] += addfac
        # TODO Still not sure if greenfac correction should be considered
        # on death side too
        curr.w *= 0.9
        curr.w[curr.k+1] = 0.1
        curr.k += 1
        curr.mu[:, curr.k+1:end] = 0.0

    elseif scheme == 1 # Death
        curr.w[curr.k] = 0.0
        corrfac = sum(curr.w[1:curr.k-1]) # TODO check -1 factor here
        curr.w /= corrfac
        curr.k -= 1
        #curr.mu[:, curr.k] = 0.5*curr.mu[:,curr.k] + 0.5*curr.mu[:, curr.k+1]
        curr.mu[:, curr.k+1:end] = 0.0
        obj.green_factor = -curr.k * log(corrfac) # Negative is because we want 1/corrfac
        
    elseif scheme == 2 # Modify params
        if curr.k != 1
            i = rand(1:curr.k)
            curr.w[1:curr.k] = curr.w[1:curr.k] + randn(curr.k) * propmoves.wmove
            curr.w = clamp(curr.w, 0.05, Inf)
            curr.w[1:curr.k] = curr.w[1:curr.k] / sum(curr.w[1:curr.k])
        else
            i=1
        end

        #Modify means
        curr.mu[:,i] += randn(prior.D) * propmoves.mumove

        #Modify di's
        #curr.d = clamp(curr.d + randn(self.n) * 0.2, 8,12) #FIXME 

        #Modify covariances
        curr.sigma = rand(InverseWishart(propmoves.priorkappa, 
                        curr.sigma*(propmoves.priorkappa-1-prior.D)))

        for i in 1:length(curr.lam) 
            curr.lam[i] += randn() * propmoves.lammove
        end
    end
    #if abs(obj.green_factor) > 2
        #println("green: $(obj.green_factor), scheme: $scheme")
    #end
    return scheme
end

import Base.sum
sum(x::Float64) = x

function energy(obj::MPMCls)
    accum = 0.0
    #likelihoods
    for i in 1:size(obj.data,1)
        for j in 1:size(obj.data,2)
            accum -= logpdf(Poisson(obj.d[i] * exp(obj.curr.lam[j,i])), 
                        floor(obj.data[i,j])) #FIXME lam and data dimensions
        end
    end
    #priors
    #lambdas
    lams = [sum(logpdf(MultivariateNormal(obj.curr.mu[:,i], obj.curr.sigma), obj.curr.lam)) 
        for i in 1:obj.curr.k]
    accum -= logsum(lams, obj.curr.w[1:obj.curr.k])

    if obj.usepriors
        # Sigma
        accum -= logpdf(InverseWishart(obj.prior.kappa, obj.prior.S), obj.curr.sigma)
        # k
        accum -= logpdf(Geometric(obj.prior.mix_comps), obj.curr.k)
        # mu
        for k in 1:obj.curr.k
            accum -= sum(logpdf(
                DiagNormal(obj.prior.mu[:,k], obj.prior.mu_sigmas), obj.curr.mu[:,k]))
        end
        # TODO w
    end
    obj.curr.energy = accum
    return accum - obj.green_factor
end
    
reject(obj::MPMCls) = copy!(obj.curr, obj.old)

function gen_points(n, d, pt)
    newpts = Array(Float64, size(pt.lam,1), n)
    # FIXME assuming single mixture
    for i=1:n
        newpts[:,i] = map(x->rand(Poisson(d*exp(x))), rand(MultivariateNormal(pt.mu[:,1], pt.sigma)))
    end

    # Or using lams below is another way to go, but it is not ideal for getting a true
    # 'feel' for the posterior fit
    #
    #for (i,lamind)=enumerate(sample(1:size(pt.lam,2), n))
        #newpts[:,i] = map(x->rand(Poisson(d*exp(x))), pt.lam[:,lamind])
    #end

    return newpts
end

function gen_posterior_points(n, d, posterior)
    npost = length(posterior)
    nper = iceil(n/npost)
    pts = Array(Float64, size(posterior[1].lam, 1), nper*npost)
    for i in 1:npost
        pts[:, (i-1)*nper+1:i*nper] = gen_points(nper, d, posterior[i])
    end
    return pts
end

function calc_pvals(Ts, points, db, d=10.0)
    numpoints = size(points, 2) # FIXME Which direction are we going again?
    ndims = size(points,1)
    tvals = Array(Float64, ndims, length(Ts), length(db))
    tru_ts = Array(Float64, ndims, length(Ts))
    pvals = similar(tru_ts)
    for i in 1:length(db)
        newpts = gen_points(numpoints, d, db[i])
        for dim=1:ndims
            for (j,T) in enumerate(Ts)
                tvals[dim,j,i] = T(vec(newpts[dim,:]))
            end
        end
    end
    for dim=1:ndims
        for (j,T) in enumerate(Ts)
            tru_ts[dim,j] = T(vec(points[dim,:]))
            pvals[dim,j] = mean(tru_ts[dim,j] .< tvals[dim,j,:]) + 0.5*mean(tru_ts[dim,j] .== tvals[dim,j,:])
        end
    end
    return pvals, tru_ts, tvals
end

function bee_moments(points,db1::Vector{Any},db2::Vector{Any},numlam,volume;dmean=10.0)
    # FIXME This is only valid for c=0.5
    g1 = calc_g(points, db1, numlam)
    g2 = calc_g(points, db2, numlam)
    bee_moments(points,db1,db2,g1,g2,numlam,volume,dmean=dmean)
end

function bee_moments(points,db1::Vector{Any}, db2::Vector{Any}, 
        g1::Vector{Float64},g2::Vector{Float64},numlam,volume;dmean=10.0)
    # FIXME This is only valid for c=0.5
    bee = exp(logsum(min(g1,g2) .+ log(volume)))*0.5
    numpts = size(points, 1)
    dims = size(points,2)
    dblen = length(db1)
    assert(length(db1) == length(db2))

    lams1 = Array(Float64, dims, numlam)
    lams2 = Array(Float64, dims, numlam)
    llams1 = Array(Float64, dims, numlam)
    llams2 = Array(Float64, dims, numlam)
    rlams1 = Array(Float64, dims, numlam)
    rlams2 = Array(Float64, dims, numlam)
    bee2 = dbpass = 0.0
    points = points'
    lgpoints = lgamma(points .+ 1)

    for i in 1:dblen 
        dbpass = 0.0
        curr1 = db1[i]
        curr2 = db2[i]
        assert(curr1.k==1 && curr2.k==1)
        rand!(MultivariateNormal(curr1.mu[:,1], curr1.sigma), lams1)
        rand!(MultivariateNormal(curr2.mu[:,1], curr2.sigma), lams2)
        for k=1:numlam*dims 
            rlams1[k] = dmean * exp(lams1[k])
            llams1[k] = log(dmean) + lams1[k]
            rlams2[k] = dmean * exp(lams2[k])
            llams2[k] = log(dmean) + lams2[k]
        end
        for j in 1:numpts # each point
            accumlam1 = 0.0
            accumlam2 = 0.0
            for s in 1:numlam # each lambda
                accumD1 = 0.0
                accumD2 = 0.0
                for d in 1:dims # each dimension
                    accumD1 += points[d,j]*llams1[d,s] - rlams1[d,s] - lgpoints[d,j]
                    accumD2 += points[d,j]*llams2[d,s] - rlams2[d,s] - lgpoints[d,j]
                end
                accumlam1 += exp(accumD1)
                accumlam2 += exp(accumD2)
            end
            temp = g1[j] < g2[j] ? accumlam1 : accumlam2
            #temp = accumlam1 # for testing purposes
            dbpass += temp
        end
        temp = 0.5*dbpass/numlam*volume
        bee2 += temp^2
    end
    bee2 /= dblen
    return bee, bee2
end

function bee_e_eff(g0, g1, volume)
    # FIXME This is only valid for c=0.5
    # Calculate expectation of error given the two effective class conditional
    # densities (posterior predictive densities)
    exp(logsum(min(g0,g1) .+ log(volume)))*0.5
end

function bee_e_grid(points, db1, db2, numlam, volume; dmean=10.0)
    # FIXME This is only valid for c=0.5
    g1 = calc_g(points, db1, numlam, dmean=dmean)
    g2 = calc_g(points, db2, numlam, dmean=dmean)
    bee_e_eff(g1,g2,volume)
end

function bee_e_data(data, db1, db2, numlam; dmean=10.0)
    # FIXME This is only valid for c=0.5
    local volume, points, g1, g2
    g1sum = g2sum = 0.0
    factor = 0.0
    D = size(data,2)
    maxN = iround(20_000^(1/D))
    trycount = 1
    while max(abs(g1sum-1.0),abs(g2sum-1.0)) > 0.1
        mins, maxs = get_bbox(data, factor=factor)
        lens, steps, points = gen_grid(mins, maxs, maxN)
        volume = prod(steps)
        g1 = calc_g(points, db1, numlam, dmean=dmean)
        g2 = calc_g(points, db2, numlam, dmean=dmean)
        g1sum = exp(logsum(g1 .+ log(volume)))
        g2sum = exp(logsum(g2 .+ log(volume)))
        println("g1 sum: $g1sum")
        println("g2 sum: $g2sum")
        #println("steps: $steps")
        #println(maxs, " ", size(points), " ", factor)
        factor += 1.0
        if trycount > 10
            return -min(g1sum, g2sum)
        else
            trycount += 1
        end
    end
    bee_e_eff(g1,g2,volume)
end

function bee_e_cube(data, db1, db2, numlam; dmean=10.0, abstol=0.03, maxevals=0)
    # FIXME This is only valid for c=0.5
    # Get the first moments using histories and adaptive integration by
    # calling Cubature.hquadrature_v with a function that evaluates the error
    # at multiple points at a time
    mins, maxs = get_bbox(data, factor=2)
    assert(length(db1) == length(db2))
    global iters = 0
    global evals = 0
    function error_1st(points, vals)
        #points: dxn array to evaluate at
        #vals: 2xn values to store into
        numpts = size(points, 2)
        iters += 1
        evals += numpts
        g1 = calc_g(points', db1, numlam, dmean=dmean)
        g2 = calc_g(points', db2, numlam, dmean=dmean)
        vals[:] = exp(min(g1,g2))
    end
    val, err = hcubature_v(error_1st, mins, maxs, abstol=abstol, maxevals=maxevals)
    println("maxs: $maxs, mins: $mins")
    println("Cubature used $iters iterations, and $evals * numclasses evaluations")
    return val*0.5, err
end

function calc_g(points, db, numlam; dmean=10.0)
    numpts = size(points, 1)
    dims = size(points, 2)
    dblen = length(db)
    res = zeros(numpts)
    lams = Array(Float64, dims, numlam)
    llams = Array(Float64, dims, numlam)
    rlams = Array(Float64, dims, numlam)
    points = floor(points)'
    lgpoints = lgamma(points .+ 1)
    for i in 1:dblen # each draw of theta
        curr = db[i]
        assert(curr.k==1)
        rand!(MultivariateNormal(curr.mu[:,1], curr.sigma), lams)
        for k=1:numlam*dims 
            rlams[k] = dmean * exp(lams[k])
            llams[k] = log(dmean)+lams[k]
        end
        for j in 1:numpts # each point
            accumlam = 0.0
            for s in 1:numlam # each lambda
                accumD = 0.0
                for d in 1:dims # each dimension
                    accumD += points[d,j]*llams[d,s] - rlams[d,s] - lgpoints[d,j]
                end
                accumlam += exp(accumD)
            end
            res[j] += (accumlam/numlam)
        end
    end
    return log(res / length(db))
end

function predict(db0, db1, points)
    g0 = calc_g(points, db0, 20)
    g1 = calc_g(points, db1, 20)
    return vec((g0 .- g1) .< 0) * 1
end

function error_points(db0, db1, points, labels)
    return sum(abs(predict(db0, db1, points) - labels))/size(labels,1)
end

end
