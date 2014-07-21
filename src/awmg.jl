type AMWGRecord <: MCMC
    obj :: Sampler
    blocks :: Int
    batchsize :: Int
    sigmas :: Vector{Float64}
    target :: Float64
    mapvalue :: Sampler
    mapenergy :: Float64
    db :: Vector{Any}

    block_accept :: Vector{Int}
    iterations :: Int

    burn :: Int
    thin :: Int
end

AMWGRecord(obj::Sampler, blocks; burn=0, thin=1) = AMWGRecord(deepcopy(obj), 
                                blocks, #Gibbs updates
                                50, # batchsize
                                ones(blocks), #sigmas
                                0.44, #target
                                deepcopy(obj), #mapvalue
                                Inf, #mapenergy
                                Any[], #db
                                zeros(Int,blocks),1, #accept, iterations
                                burn,thin) #burn, thin

import Distributions.sample

function sample(rec::AMWGRecord, iters::Int; verbose=false, adjust="burnin")
    oldenergy = energy(rec.obj)
    if verbose
        println("Initial energy: $oldenergy")
    end

    for current_iter = rec.iterations:(rec.iterations+iters-1)
        sigma = randn(rec.blocks) .* rec.sigmas
        for i in 1:rec.blocks
            save!(rec.obj)
            #propose new object in gibbs block i, with sigma sigma[i]
            propose!(rec.obj, i, rec.sigmas[i])
            newenergy = energy(rec.obj, i) 
            if newenergy < rec.mapenergy
                rec.mapvalue = deepcopy(rec.obj)
                rec.mapenergy = newenergy
            end
            r = oldenergy - newenergy
            ### Acceptance of new moves ###
            if rand() < exp(r)
                rec.block_accept[i] += 1
                oldenergy = newenergy
            else
                reject!(rec.obj)
            end
        end

        # Adjust weights
        if current_iter % rec.batchsize == 0 && ((adjust == "burnin" && current_iter <= rec.burn) || adjust == "always")
            delta = min(0.1, (current_iter / rec.batchsize)^-0.5)
            verbose && println("Sigmas before tuning: $(rec.sigmas)")
            verbose && println("Delta: $delta")
            verbose && println("Block_accept: $(rec.block_accept)")
            for i in 1:length(rec.sigmas)
                rec.sigmas[i] *= rec.block_accept[i] / rec.batchsize < rec.target ? exp(-delta) : exp(delta)
                rec.block_accept[i] = 0 # FIXME, is this right? I couldn't find it in mamba, but ow above doesn't make sense
            end
            verbose && println("Sigmas after tuning: $(rec.sigmas)")
        end

        if rec.iterations >= rec.burn && rec.iterations%rec.thin == 0
            push!(rec.db, deepcopy(rec.obj.curr))
        end
        
        if (rec.iterations) % 1000 == 0 && verbose
            @printf "Iteration: %8d, best energy: %7f, current energy: %7f\n" rec.iterations rec.mapenergy oldenergy
        end
        rec.iterations += 1
    end
    if verbose
        println("Accepted samples: $(rec.block_accept)")
        println("Acceptance ratios: $(rec.block_accept./(rec.iterations-rec.burn))")
    end
end

# To be overwritten
energy(x::Sampler) = error("energy() must be instantiated for your sampler object")
propose!(x::Sampler) = error("propose!() must be instantiated for your sampler object")
reject!(x::Sampler) = error("reject!() must be instantiated for your sampler object")
save!(x::Sampler) = error("save!() must be instantiated for your sampler object")
