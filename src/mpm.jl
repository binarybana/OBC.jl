module MPM
    
using Distributions
using OBC
using SAMC
import SAMC: energy, propose!, save!, reject!, record

export MPMPrior, MPMParams, MPMSampler, calc_g, gen_points, 
    gen_posterior_points, calc_pvals, error_points, predict, error_moments_cube, 
    var_error_hists, e_error_hists, e_error_eff, mpm_classifier, sample!

include("mpmcls.jl")

######################################################################
######################################################################
######################################################################
######################################################################

type MPMPrior
    D :: Int
    mu :: Vector{Float64}
    mu_sigmas :: Vector{Float64} # Unvariate variances on mu
    kappa :: Float64
    S :: Matrix{Float64}
end

function MPMPrior(;D=2, args...) 
    kappa = 6.0
    p = MPMPrior(D, zeros(D), ones(D)*30, kappa, (kappa-1-D)*eye(D,D))
    for (sym, val) in args
        if sym in fieldnames(MPMPrior)
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
    mu :: Vector{Float64}
    sigma :: Matrix{Float64}
    sigpre :: Matrix{Float64}
    lam :: Matrix{Float64}
    energy :: Float64 # I dunno... does this belong here?
    # d :: Vector{Float64}
    function MPMParams(mu, sigpre, lam)
        @assert size(sigpre) == (size(mu,1),size(mu,1))
        @assert size(lam,1) == size(mu,1)
        new(mu, sigpre'*sigpre, sigpre, lam, -Inf)
    end
end

function copy!(to::MPMParams, from::MPMParams)
    to.mu[:] = from.mu
    to.sigma[:] = from.sigma
    to.sigpre[:] = from.sigpre
    to.lam[:] = from.lam
    to.energy = from.energy
end

######################################################################
######################################################################
######################################################################
######################################################################

type MPMSampler <: Sampler
    curr :: MPMParams
    old :: MPMParams
    prior :: MPMPrior
    data :: Matrix{Int}
    usepriors :: Bool
    d :: Vector{Float64}
end

function MPMSampler(prior::MPMPrior, data::Matrix, obj::MPMParams, d::Float64)
    MPMSampler(prior, data, obj, d*ones(size(data,1)))
end

function MPMSampler(prior::MPMPrior, data::Matrix, obj::MPMParams, d::Vector{Float64})
    @assert size(data,1) > size(data,2)
    if !(eltype(data) <: Integer)
        warn("Data is being truncated from type $(eltype(data)) to Int, dataloss might occur")
    end
    MPMSampler(deepcopy(obj), deepcopy(obj), prior, round(Int, data), true, d)
end

const blocks = 3

function propose!(obj::MPMSampler, block::Int, sigma::Float64)
    curr = obj.curr
    prior = obj.prior

    if block == 1 #Modify means
        curr.mu[:] += randn(prior.D) * sigma*0.1

    elseif block == 2

        #Modify covariances
        #temp = curr.sigma .* (0.5 .+ rand(prior.D, prior.D)) #.* curr.sigma ./ mean(curr.sigma)
        #temp = curr.sigma .+ diagm(abs(clamp(diag(curr.sigma),0.001,Inf)) .* (rand(prior.D).-0.5))
        #temp = (temp'*temp).^(1/2)
        #if isposdef(temp) && false #rand(1:100) > 2
            #curr.sigma = temp
        #else
        for i=1:prior.D*prior.D
            curr.sigpre[i] += rand(Normal(0.0,sigma*0.1))
        end
        At_mul_B!(curr.sigma, curr.sigpre, curr.sigpre) # sigma = sigpre'*sigpre

        #curr.sigma = rand(InverseWishart(propmoves.priorkappa, 
                #curr.sigma*(propmoves.priorkappa-1-prior.D)))
        #end

    elseif block == 3
        # Modify lambdas
        for i in 1:length(curr.lam) 
            if rand(Bool)
                curr.lam[i] += randn() * sigma*0.1
            end
        end

    #elseif block == 3
        #Modify di's
        #curr.d = clamp(curr.d + randn(self.n) * 0.2, 8,12) #FIXME 

    end
    nothing
end

function energy(obj::MPMSampler, block::Int=0) #block currently unused
    accum = 0.0
    #likelihoods
    for i in 1:size(obj.data,1), j in 1:size(obj.data,2)
        accum -= logpdf(Poisson(obj.d[i] * exp(obj.curr.lam[j,i])), 
                    obj.data[i,j])
    end
    #lambdas
    accum -= sum(logpdf(MultivariateNormal(obj.curr.mu, obj.curr.sigma), obj.curr.lam)) 

    if obj.usepriors
        # Sigma
        #accum -= sum(logpdf(Normal(0.0, obj.propmoves.premove), obj.sigpre)) 
        accum -= logpdf(InverseWishart(obj.prior.kappa, obj.prior.S), obj.curr.sigma)
        # mu
        accum -= sum(logpdf(
            MvNormal(obj.prior.mu, obj.prior.mu_sigmas), obj.curr.mu))
    end
    obj.curr.energy = accum
    return accum
end

record(obj::MPMSampler) = deepcopy(obj.curr)
    
reject!(obj::MPMSampler) = copy!(obj.curr, obj.old)
save!(obj::MPMSampler) = copy!(obj.old, obj.curr)

function gen_points(n, dmean, pt)
    D = size(pt.lam,1)
    newpts = Array(Float64,D,n)
    for i=1:n
        tempmv = rand(MultivariateNormal(pt.mu, pt.sigma))
        for j=1:D
            rate = dmean*exp(tempmv[j])
            if rate > 1000
                newpts[j,i] = rand(Normal(rate,sqrt(rate)))
            else
                try
                    newpts[j,i] = rand(Poisson(rate))
                catch x
                    @show newpts[j,i]
                    @show rate
                    throw(x)
                end
            end
        end
    end

    # Or using lams below is another way to go, but it is not ideal for getting a true
    # 'feel' for the posterior fit
    #
    #for (i,lamind)=enumerate(sample(1:size(pt.lam,2), n))
        #newpts[:,i] = map(x->rand(Poisson(d*exp(x))), pt.lam[:,lamind])
    #end

    return newpts
end

function gen_posterior_points(n, dmean, posterior)
    npost = length(posterior)
    nper = iceil(n/npost)
    pts = Array(Float64, size(posterior[1].lam, 1), nper*npost)
    for i in 1:npost
        pts[:, (i-1)*nper+1:i*nper] = gen_points(nper, dmean, posterior[i])
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

function calc_g(points, db, numlam; dmean=10.0)
    # Note: dmean here could actually be d if we knew valid values of d for
    # each point we are calculating the effective density for
    numpts = size(points, 2)
    dims = size(points, 1)
    dblen = length(db)
    res = zeros(numpts)
    lams = Array(Float64, dims, numlam)
    llams = Array(Float64, dims, numlam)
    rlams = Array(Float64, dims, numlam)
    points = floor(points)
    lgpoints = lgamma(points .+ 1)
    for i in 1:dblen # each draw of theta
        curr = db[i]
        rand!(MultivariateNormal(curr.mu, curr.sigma), lams)
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

# Utility convenience function
function mpm_model(data; burn=1000, thin=50, d=100.0, usepriors=true)
    D = size(data, 2)
    cov = eye(D,D) .* 0.1
    mu = ones(D)
    prior = MPM.MPMPrior(D=D, kappa=10., S=eye(D)*0.05*10)

    start = MPM.MPMParams(mu, cov, 
        clamp(log(data'/d),-8.0,Inf)) #lam 
    obj = MPM.MPMSampler(prior, data, deepcopy(start), d)
    obj.usepriors = usepriors
    #mymh_a = MPM.MHRecord(obj_a,burn=burn,thin=thin)
    return MPM.AMWGRecord(obj, blocks, burn=burn,thin=thin)
end

end
