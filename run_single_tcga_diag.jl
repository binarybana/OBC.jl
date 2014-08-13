reload("src/OBC.jl")
using HDF5, JLD
using Stats
using Distributions

loc = expanduser("~/largeresearch/seq-data/rnaseq.jld")
@load loc luad_norm luad_norm_sampids luad_norm_geneids

blas_set_num_threads(1)

d = 10.0
iters = 1000
burn = 200
thin = 50
D = 2
rounds = 1
numsamps = 100

samps = sample(1:size(luad_norm,2), numsamps)
genes = sample(1:size(luad_norm,1), D*rounds)

i=1
data = luad_norm[genes[(i-1)*D+1:i*D],samps]'

######################################################################
# Fit Model
######################################################################

#samplers = MPM.AMWGRecord[]
#for i=1:5
    #push!(samplers, MPM.mpm_model(data,burn=burn,thin=thin,d=d,usepriors=true))
    #sample(samplers[end],iters)
#end

import Base: length
length(x::MPM.AMWGRecord) = length(x.db)

function gelman_rubin(samplers)
    # parameters in db must be scalars, vectors or matrices
    assert(all((map(length,samplers) .- length(samplers[1])).==0)) # ... for now

    m = length(samplers)
    n = length(samplers[1])

    paramnames = names(samplers[1].db[1])
    dnostics = Dict()
    for p in paramnames
        ex = getfield(samplers[1].db[1],p)
        extype = typeof(ex)


        if extype <: Number
            posts = [convert(Vector{extype}, map(y->getfield(y,p), x.db)) for x in samplers]
        elseif ndims(ex) == 1
            posts = vcat([map(y->getfield(y,p), x.db) for x in samplers]...)
        elseif ndims(ex) == 2
            posts = vcat([map(y->vec(getfield(y,p)), x.db) for x in samplers]...)
        end

        vars = vcat(map(x->var(x,1), posts)...)
        W = mean(vars,1)
        overallmean = mean(map(mean,posts))
        B = n/(m-1).*sum((vars .- overallmean).^2)
        dnostics[p] = sqrt(((n-1)/n .* W .+ B./n) ./ W)
    end
    return dnostics
end

dnostics = gelman_rubin(samplers)
