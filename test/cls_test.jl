
module ClsTest

using Distributions
using Sundials
using PyPlot

abstract Sampler

include("mh.jl")

SAMPS = 10
include("cls_utils.jl")

type DynPrior
    theta_mu :: Float64
    theta_sigma :: Float64
    init1_mu :: Float64
    init1_sigma :: Float64
    init2_mu :: Float64
    init2_sigma :: Float64
    sigma_alpha :: Float64
    sigma_beta :: Float64
end
DynPrior(theta_mu::Float64) = DynPrior(theta_mu, 
    10.0, #theta_sigma :: Float64
    0.0, #init1_mu :: Float64
    10.0, #init1_sigma :: Float64
    1.0, #init2_mu :: Float64
    10.0, #init2_sigma :: Float64
    2.0, #sigma_alpha :: Float64
    1.0) #sigma_beta :: Float64

type DynParams 
    theta :: Float64
    sigma :: Float64
    init1 :: Float64
    init2 :: Float64
end

draw_params(p::DynPrior) = DynParams(rand(Normal(p.theta_mu, p.theta_sigma)),
                                    rand(InverseGamma(p.sigma_alpha, p.sigma_beta)),
                                    rand(Normal(p.init1_mu, p.init1_sigma)),
                                    rand(Normal(p.init2_mu, p.init2_sigma)))

function draw_data(params::DynParams, n=10)
    currpred = pred(params.theta, params.init1, params.init2)
    preds = repmat(currpred', n)
    preds += rand(Normal(0.0, params.sigma), size(preds)...)
    #preds += rand(Normal(0.0, params.sigma), size(preds)...)
    #preds = preds .+ linspace(0,0.5,SAMPS)'
end

type DynCls <: Sampler
    curr :: DynParams
    old :: DynParams
    prior :: DynPrior
    data :: Matrix{Float64}
end

DynCls(prior::DynPrior, data::Matrix{Float64}) = DynCls(draw_params(prior), draw_params(prior), prior, data)
DynCls(prior::DynPrior, data::Matrix{Float64}, obj::DynParams) = DynCls(obj, deepcopy(obj), prior, data)

function propose(obj::DynCls)
    copy!(obj.old, obj.curr)
    obj.curr.theta += randn() * 0.1
    obj.curr.sigma += randn() * 0.08
    obj.curr.sigma = clamp(obj.curr.sigma, 0.01, Inf)
    obj.curr.init1 += randn() * 0.1
    obj.curr.init2 += randn() * 0.1
    0
end

function pred(theta, init1, init2)
    #Have to use a closure because the API is missing user data atm:
    function f(t,y,ydot)
        ydot[1] = y[2]
        ydot[2] = -theta*y[1]
    end
    Sundials.ode(f, [init1, init2], linspace(0,2,SAMPS))[:,1]
end

function energy(obj::DynCls)
    accum = 0.0
    currpred = pred(obj.curr.theta, obj.curr.init1, obj.curr.init2)
    #likelihoods
    accum -= sum(logpdf(Cauchy(0.0, obj.curr.sigma), vec(obj.data .- currpred')))
    #accum -= sum(vec(obj.data .- currpred').^2)
    #accum -= sum(logpdf(Normal(0.0, obj.curr.sigma), vec(obj.data .- currpred')))
    #priors
    accum -= logpdf(Normal(obj.prior.theta_mu, obj.prior.theta_sigma), obj.curr.theta) #theta
    accum -= logpdf(InverseGamma(obj.prior.sigma_alpha, obj.prior.sigma_beta), obj.curr.sigma) #sigma
    accum -= logpdf(Normal(obj.prior.init1_mu, obj.prior.init1_sigma), obj.curr.init1) #init1
    accum -= logpdf(Normal(obj.prior.init2_mu, obj.prior.init2_sigma), obj.curr.init2) #init2
end
    
reject(obj::DynCls) = copy!(obj.curr, obj.old)

function copy!(to::DynParams, from::DynParams)
    to.theta = from.theta
    to.sigma = from.sigma
    to.init1 = from.init1
    to.init2 = from.init2
end

function mc_approx(traces, db)
    accum = zeros(size(traces,1), length(db))
    for (i,params) in enumerate(db)
        currpred = pred(params.theta, params.init1, params.init2)
        for j in 1:size(traces,1)
            accum[j,i] = exp(sum(logpdf(Cauchy(0.0, params.sigma), traces[j,:]' .- currpred)))
            #accum[j,i] = sum(vec(traces[j,:]' .- currpred).^2)
        end
    end
    mapslices(sum, accum, 2)
end

function classify(traces::Matrix{Float64}, db_a, db_b; labels=None)
    pxy = [mc_approx(traces, db_a) mc_approx(traces, db_b)]
    preds = mapslices(x->findmax(x)[2], pxy, 2)
    logdiffs = log(pxy[:,1]) - log(pxy[:,2])
    err = 0
    for i in 1:length(preds)
        if labels != None
            #println("Cls prediction: $(preds[i]), true label: $(labels[i]), log difference: $(logdiffs[i])")
            if preds[i] != labels[i]
                err += 1
            end
        else
            #println("Cls prediction: $(preds[i]), log difference: $(logdiffs[i])")
        end
    end
    if labels != None
        println("Error rate: $(err/length(preds))")
    end
    preds
end

prior_a = DynPrior(40.0)
prior_b = DynPrior(30.0)

start_both = DynParams(30.0, 0.5, 1.0, 0.9)
true_a = DynParams(35.0, 0.5, 1.3, 1.2)
true_b = DynParams(30.0, 0.75, 1.0, 1.0)
#true_a = draw_params(prior_a)
#true_b = draw_params(prior_b)

data_a = draw_data(true_a)
data_b = draw_data(true_b)

samp_a = DynCls(prior_a, data_a, start_both)
samp_b = DynCls(prior_b, data_b, start_both)

mymh_a = MHRecord(samp_a,burn=2000,thin=20)
mymh_b = MHRecord(samp_b,burn=2000,thin=20)

sample(mymh_a,5000)
sample(mymh_b,5000)

plot_all(mymh_a,mymh_b,hidden=true)
legend()
plot(samp_a.data[1,:]',"k--", linewidth=3)
plot(samp_b.data[1,:]',"k.-", linewidth=3)

#plot_all(mymh_a,mymh_b)
#legend()
#plot(samp_a.data[1,:]',"k--")
#plot(samp_b.data[1,:]',"k.-")

classify([draw_data(true_a,1000), draw_data(true_b,1000)], mymh_a.db, mymh_b.db, labels=[ones(Int,1000), ones(Int,1000)+1])

#plot_trace(mymh_a, mymh_b, "sigma")
#plot_trace(mymh_a, mymh_b, "theta")
#plot_posterior(mymh_a)
#plot_posterior(mymh_b)
#
# Talk about:
# Model formulation: x'' = -\theta x x(0) = init1, x'(0)=init2
# Cauchy vs Normal
# convergence and starting point
# epsilon tube
#

end
